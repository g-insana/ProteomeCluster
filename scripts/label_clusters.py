#!/usr/bin/env python3
# changelog
# Wed  6 Nov 2024 10:05:27 GMT 0.1 started
# Wed  6 Nov 2024 10:05:27 GMT 1.0 working sequential approach
# Wed  6 Nov 2024 20:15:49 GMT 1.1 added --nolabel
# Wed  6 Nov 2024 21:44:22 GMT 1.8 coded P1 of parallel approach
# Thu  7 Nov 2024 12:32:01 GMT 1.9 coded P2
# Thu  7 Nov 2024 15:55:41 GMT 2.0 working parallel approach
# Sun 10 Nov 2024 16:26:30 GMT 2.1 alternative parallel approach (nochunks)
# Sun 10 Nov 2024 18:30:42 GMT 2.2 rewrote docstrings
# Tue 12 Nov 2024 17:39:07 GMT 2.3 added timeout and retries for file locking
# Wed 13 Nov 2024 11:56:31 GMT 2.4 added possibility to work in batches but single-thread
# Sat 30 Nov 2024 16:02:39 GMT 2.5 added progress bar
# Sun 23 Nov 17:54:06 GMT 2025 2.6 compute_split_positions and split_file from ffdb

# imports
import os
import sys
import re
import time
from tqdm import tqdm
from ffdb import compute_split_positions, split_file
import shutil
import argparse
from glob import glob
from typing import List, Optional


# for distributed/batch approach only:
from random import sample, randint
from multiprocessing import Pool, current_process, cpu_count
from filelock import Timeout, FileLock
from math import ceil

# for profiling/testing only:
# from pympler import asizeof

# constants
READMETHOD = "lines"
HEADER = "cluster_id\tprotein_id\tproteomes\tis_rep\n"
BUFFERSIZE = 1048576  # 1Mb
MINCHUNKSIZE = "5m"  # 5 Mb
RE_REMOVE_EXTENSION = re.compile(r"\.[^.]+$")
DESCRIPTION = """
This script labels the second column of a TSV (protein IDs) with the FASTA filename(s) (proteome identifiers) where each protein appears, compacts cluster IDs to sequential integers, and flags the representative protein for each cluster. Input row order is preserved.

 * Expects a two-column TSV (like mmseqs2 cluster files): column 1 = cluster identifier, column 2 = protein identifier (matching FASTA headers).
 * Scans FASTA files in FASTA_DIR and maps protein IDs (from headers) to the FASTA filename
 * Adds a new column, proteomes, containing the proteome identifier(s) for each protein. If a protein is found in multiple FASTA files, identifiers are comma-separated.
 * Replaces original cluster identifiers with sequential integers starting at 0 while preserving the original row order and grouping identical original cluster IDs to the same integer.
 * Adds a fourth column, is_rep, marking the first protein encountered for each original cluster with '*' (representative); other rows in that cluster have an empty value.
 * Produces an output TSV with header:
    cluster_id\tprotein_id\tproteomes\tis_rep

Required arguments:
    FASTA_DIR  - directory containing FASTA files to search (headers provide protein IDs)
    INPUT_FILE - input TSV file (two columns: cluster identifier and protein identifier)
    OUT_FILE   - output TSV to create (four columns)

Sample input file (cluster and protein columns (tab-separated):
    ENSSSCP00000055324|661\tENSSSCP00000055324|661
    ENSSSCP00000055324|661\tENSSSCP00055011301|568
    ENSSSCP00000055324|661\tENSSSCP00035021320|596

Sample output file (header included):
    cluster_id\tprotein_id\tproteomes\tis_rep
    0\tENSSSCP00000055324|661\t35497\t*
    0\tENSSSCP00055011301|568\t4698922\t
    0\tENSSSCP00035021320|596\t4698918\t

Example call:
    label_clusters.py --fasta_dir pig/ --input_file results_pig/Specie_protein_cluster.tsv --out_file results_pig/Labelled_Specie_protein_cluster.tsv --prefix proteome_ --extension .fa

Distributed approach (batch tagging, optionally parallel):
      When dealing with a huge amount of entries memory requirements could be a problem.
      The labelling work could then be performed in batches, specifying a --batchsize. That will be the
      maximum number of fasta_files to read in one go to create a partial dictionary and annotate the
      input file with that. The input file will be annotated several times until all fasta_files have
      been ingested.

      This can be combined with --threads to do the tagging in parallel.
      The input file will be split into thrice the number of threads and each thread will only read
      a number of fasta files equal to batchsize at each time.

      Optionally --chunksize can be specified: the input file will be split in chunks of that size
          (minimum 5Mb)

      E.g.: ./label_clusters.py --fasta_dir ecoli/ --input_file results_ecoli/Specie_protein_cluster.tsv --out_file results_ecoli/Labelled_Specie_protein_cluster.tsv --prefix proteome_ --extension .fa --threads 10 --batchsize 1000 --chunksize 1G

      Alternative approach: the input file will be copied and processed independently by each worker and
      results combined in the end. To use this approach specify --chunksize=n.
      This approach eliminates the chance of time spent with workers trying to acquire file lock on
      the same chunk, but it uses more disk space and could suffer i/o performance loss when combining
      result files in the end.
"""

# timing tests
#    on a node in a linux datacentre:
#        1s to process 13 files of approx 30M each for a total of 627957 ids
#          'chunks': ~800000 ids/s
#          'full': ~770000 ids/s
#          'lines': ~1100000 ids/s
#        (on macos 'full' was as fast as 'lines' and 'chunks' slightly less)
#        1h30m to process ~120k files of approx 1.4Mb each for a total of 591 million ids
#          'full': ~125000 ids/s
#          'lines': ~110000 ids/s
#
# memory tests
#        55.82GiB to process ~120k files of approx 1.4Mb each for a total of 591 million ids sequentially


# helper functions
def secs2time(secs):
    """
    Converts a time duration in seconds to a human-readable format (hours, minutes, seconds).

    Args:
        secs (int): Time duration in seconds.

    Returns:
        str: A formatted string representing the time in "HHh MMm SSs" format.

    Example:
        secs2time(3663)  # Output: "01h 01m 03s"
    """
    minutes, seconds = divmod(secs, 60)
    hours, minutes = divmod(minutes, 60)
    return "{:02.0f}h {:02.0f}m {:02.0f}s".format(hours, minutes, seconds)


def elapsed_time(start_time, work_done=None):
    """
    Computes the elapsed time from a given start time in seconds and returns a formatted string.
    If `work_done` is specified, also computes the speed of the process.

    Args:
        start_time (float): The start time in seconds (from time.time()).
        work_done (int, optional): Number of completed iterations or tasks.

    Returns:
        str or tuple: A formatted string with elapsed time if `work_done` is None,
                      otherwise a tuple with formatted elapsed time and computed speed in "it/s".

    Example:
        start_secs = time.time()
        time.sleep(2)
        print(" '-- Elapsed: {} --'".format(elapsed_time(start_secs)))
        # Output example: " '-- Elapsed: 00h 00m 02s --'"

        iterations_done = 10
        print(" '-- Elapsed: {}, {} it/s --'".format(*elapsed_time(start_secs, iterations_done)))
        # Output example: " '-- Elapsed: 00h 00m 02s, 5.0 it/s --'"
    """
    process_time = time.time() - start_time
    if work_done is None:
        return secs2time(process_time)
    process_speed = round(work_done / process_time, 2)
    return secs2time(process_time), process_speed


def exit_with_error(message: str, code: int = 1):
    """
    Prints an error message to stderr and exits the program with the specified exit code.

    Args:
        message (str): The error message to display.
        code (int): The exit code to return upon termination (default: 1).
    """
    eprint(f"   => {message}")
    sys.exit(code)


def eprint(*myargs, **kwargs):
    """
    Prints the provided arguments to stderr, useful for logging errors or status without cluttering stdout.

    Args:
        *myargs: Variable length argument list, elements to be printed.
        **kwargs: Arbitrary keyword arguments (e.g., end='\n').

    Returns:
        None
    """
    print(*myargs, file=sys.stderr, **kwargs)


def read_file_in_chunks(file_path, chunk_size=BUFFERSIZE):
    """
    Reads a file in chunks and yields only complete lines to avoid splitting lines across chunks.

    Args:
        file_path (str): Path to the file.
        chunk_size (int): Number of bytes to read per chunk (default: BUFFERSIZE).

    Yields:
        str: Full lines from the file, stripped of newline characters.

    Example:
        for line in read_file_in_chunks("file.txt", chunk_size=1024):
            print(line)
    """
    with open(file_path, "r") as fh:
        buffer = ""
        while True:
            chunk = fh.read(chunk_size)
            if not chunk:
                break  # end of file

            buffer += chunk
            lines = buffer.splitlines(True)

            # process all lines except the last one which may be incomplete
            for line in lines[:-1]:
                yield line.strip()

            # keep the remaining possibly incomplete line in the buffer
            buffer = lines[-1]

        if buffer:
            yield buffer.strip()


def delete_files(filenames: List[str], path: Optional[str] = None):
    """
    Deletes specified temporary files.

    Args:
        filenames (List[str]): List of filenames to delete.
        path (Optional[str]): Optional base directory to prepend to each filename.

    Returns:
        None
    """
    for filename in filenames:
        if path is not None:
            filename = os.path.join(path, filename)
        if os.path.isfile(filename):
            os.remove(filename)


def calculate_blocksize(filename, number_of_chunks):
    """
    Calculate the size of a single block based on the number of chunks desired and the file's size.

    Args:
        filename (str): Path to the input file.
        number_of_chunks (int): Number of chunks to divide the file into.

    Returns:
        int: The calculated size of a single block (in bytes).
    """
    filesize = os.path.getsize(filename)
    if filesize == 0:
        raise RuntimeError("File '{}' is empty".format(filename))
    blocksize = int(ceil(filesize / number_of_chunks))
    number_of_chunks = int(ceil(filesize / blocksize))
    return blocksize


def siprefix2num(numberstring):
    """
    Converts a metric prefix suffixed string (e.g., "10k", "5m") to its corresponding byte value.

    Args:
        numberstring (str): A string representing a number with a metric prefix (e.g., '10k', '1G').

    Returns:
        int: The integer value of the number represented in bytes. If no valid metric prefix is found,
             the string is converted to an integer as is.
    """
    if numberstring == "0":
        return 0

    prefix = {"k": 1024, "m": 1024**2, "g": 1024**3, "t": 1024**4, "p": 1024**5}
    numberstring = numberstring.lower()
    try:
        return int(numberstring[:-1]) * prefix[numberstring[-1]]
    except KeyError:
        # no or unknown meter-prefix
        return int(numberstring)


# functions
def check_args(DESCRIPTION):
    """
    parse arguments and check for error conditions

    Args:
        DESCRIPTION (str): Description of the program to be printed with help text.
    """

    def positive_integer(value):
        try:
            value = int(value)
            if value <= 0:
                raise argparse.ArgumentTypeError(
                    "{} is not a positive integer".format(value)
                )
        except ValueError:
            raise Exception("{} is not an integer".format(value))
        return value

    def is_valid_file(path):
        """Check if the given path is a valid, readable file."""
        if not path:
            raise argparse.ArgumentTypeError(f"File path cannot be empty or None.")

        if not os.path.exists(path):
            raise argparse.ArgumentTypeError(f"The file '{path}' does not exist.")

        if not os.path.isfile(path):
            raise argparse.ArgumentTypeError(f"The path '{path}' is not a valid file.")

        if not os.access(path, os.R_OK):
            raise argparse.ArgumentTypeError(f"The file '{path}' is not readable.")

        return path

    class CustomArgumentParser(argparse.ArgumentParser):
        def print_help(self, *args, **kwargs):
            """
            print custom text before the default help message
            """
            print(DESCRIPTION)
            super().print_help(*args, **kwargs)

    parser = CustomArgumentParser(
        description="Proteome labeller for clusters of protein identifiers."
    )
    parser.add_argument(
        "-f",
        "--fasta_dir",
        type=str,
        required=True,
        help="Directory containing all .fa files",
    )
    parser.add_argument(
        "-i",
        "--input_file",
        type=is_valid_file,
        required=True,
        help="Path to the input file",
    )
    parser.add_argument(
        "-o", "--out_file", type=str, required=True, help="Path to the output file"
    )
    parser.add_argument(
        "-p",
        "--prefix",
        type=str,
        required=False,
        help="Optionally prefix for filenames, which will be removed; e.g. 'proteome_'",
    )
    parser.add_argument(
        "-e",
        "--extension",
        type=str,
        required=False,
        help="Extension for files in fasta_dir (e.g. '.fa')",
    )
    parser.add_argument(
        "-n",
        "--nolabel",
        type=str,
        required=False,
        help="Optional string for missing mapping; default is '' (e.g. '?')",
        default="",
    )
    parser.add_argument(
        "-s",
        "--sortlabels",
        action="store_true",
        required=False,
        help="Optionally sort the labels in the output file",
        default=False,
    )
    parser.add_argument(
        "-u",
        "--uniq",
        action="store_true",
        required=False,
        help="Optionally uniq the input file (in case it has repeated lines)",
        default=False,
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=positive_integer,
        required=False,
        default=1,
        help="Number of threads for parallel processing",
    )
    parser.add_argument(
        "-b",
        "--batchsize",
        type=positive_integer,
        required=False,
        help="Max number of fasta files that will be processed at the same time",
    )
    parser.add_argument(
        "-c",
        "--chunksize",
        type=str,
        required=False,
        help=f"Chunk size in which to split the input file; if not specified: split the input file in a number of chunks equal to thrice the number of threads; minimum chunksize: {MINCHUNKSIZE}; use 'n' to avoid splitting the input file",
    )
    parser.add_argument(
        "-q",
        "--progress",
        action="store_true",
        required=False,
        help="Show a progress bar",
        default=False,
    )
    args = parser.parse_args()
    eprint(f" |-- input_file: {args.input_file}")

    # check for out_file presence and writing ability
    if os.path.isfile(args.out_file):
        exit_with_error(
            f"ERROR: File '{args.out_file}' exists; remove it to proceed.", 17
        )
    try:
        with open(args.out_file, "w"):
            pass
    except PermissionError:
        exit_with_error(f"ERROR: Cannot write to file '{args.out_file}'.", 1)
    delete_files([args.out_file])  # delete the just touched file
    eprint(f" |-- out_file: {args.out_file}")

    if not os.path.isdir(args.fasta_dir):
        exit_with_error(f"ERROR: No such directory '{args.fasta_dir}'", 2)

    # check for files in fasta_dir
    search_pattern = (
        f"{args.prefix}*{args.extension}"
        if args.extension and args.prefix
        else f"*{args.extension or ''}"
    )
    fasta_files = glob(os.path.join(args.fasta_dir, search_pattern))
    if not fasta_files:
        exit_with_error(
            f"ERROR: No matching '{args.extension or ''}' files in '{args.fasta_dir}'.",
            2,
        )

    eprint(f" |-- fasta_dir: {args.fasta_dir}")
    eprint(f" |-- fasta_files: {len(fasta_files)}")

    # validate arguments for distributed/batch approach
    args.distributed = False
    if args.batchsize:
        args.distributed = True

    # chunksize processing
    if args.chunksize and not args.distributed:
        exit_with_error(
            "ERROR: chunksize requires distributed execution. Specify --batchsize and optionally --threads",
            22,
        )
    if args.chunksize is None or args.chunksize == 0:
        # if unspecified, calculate from number of threads
        args.chunksize = calculate_blocksize(args.input_file, 3 * args.threads)
        args.nochunks = False
    else:
        args.nochunks = args.chunksize == "n"  # if 'n' then nochunks
        if not args.nochunks:
            # compute from specified number
            args.chunksize = max(
                siprefix2num(args.chunksize), siprefix2num(MINCHUNKSIZE)
            )

    if args.threads > cpu_count():
        args.threads = cpu_count()
        eprint(f" |-- WARNING: only {args.threads} threads available")
    if args.distributed:
        eprint(f" |-- threads: {args.threads}")
        eprint(f" |-- batchsize: {args.batchsize}")
        if args.nochunks:
            eprint(f" |-- chunksize: 0 (input file won't be split)")
        else:
            eprint(f" |-- chunksize: {args.chunksize} bytes")

    return args, fasta_files


def create_proteome_protein_map(fasta_files, args, read_method=READMETHOD):
    """
    Creates a dictionary mapping protein identifiers (from fasta headers)
    to proteome identifiers (from file names) by reading through specified fasta files
    in one of three methods: 'lines', 'full', or 'chunks'.

    Args:
        fasta_files (list): List of paths to the fasta files.
        args (argparse.Namespace): Parsed arguments including prefix and extension options.
        read_method (str): Method to read files; options are 'lines', 'full', and 'chunks'.

    Returns:
        dict: A dictionary mapping protein identifiers (str) to proteome identifiers (str).
    """
    proteome_protein_map = {}

    def _process_protein_id(line):
        # store protein_id extracted from fasta headers
        protein_id = line.rstrip("\n")[1:].split(" ")[
            0
        ]  # remove initial > and anything after a space
        if protein_id not in proteome_protein_map:
            proteome_protein_map[protein_id] = proteome_id
        else:
            proteome_protein_map[protein_id] += "," + proteome_id

    if read_method not in {"lines", "full", "chunks"}:
        eprint(f"ERROR: Invalid read_method '{read_method}' specified.")
        return proteome_protein_map

    for file in fasta_files:
        proteome_id = os.path.basename(file)

        # remove extension, if any
        if args.extension is not None:
            proteome_id = proteome_id[0 : -len(args.extension)]
        else:
            proteome_id = RE_REMOVE_EXTENSION.sub("", proteome_id)

        # optionally remove prefix from the filenames
        if args.prefix is not None:
            proteome_id = proteome_id[len(args.prefix) :]

        if read_method == "lines":  # read line by line
            with open(file, "r") as fh:
                for line in fh:
                    if line.startswith(">"):
                        _process_protein_id(line)
        elif read_method == "full":  # read whole file in memory
            with open(file, "r") as fh:
                for line in fh.readlines():
                    if line.startswith(">"):
                        _process_protein_id(line)
        elif read_method == "chunks":  # read file in chunks
            for line in read_file_in_chunks(file):
                if line.startswith(">"):  # fasta header
                    _process_protein_id(line)
        else:
            eprint(f"    => ERROR: no such read_method {read_method}")

    # eprint(list(proteome_protein_map.items())[:5]) #debug, first 5 items in the map
    # eprint(list(proteome_protein_map.items())[-5:]) #debug, last 5 items in the map
    # eprint(f" |-- proteome_protein_map: {len(proteome_protein_map)} keys, {asizeof.asizeof(proteome_protein_map)} bytes") #debug
    return proteome_protein_map


def label_proteins(
    input_file,
    output_file,
    proteome_protein_map,
    nolabel="",
    sortlabels=False,
    uniq=False,
):
    """
    Use the provided mapping to label the input_file, assigning proteome labels to protein_ids.
    Optionally sort (and uniq) the attached labels.

    Args:
        input_file (str): Path to the input file with protein identifiers.
        output_file (str): Path to the output file where results are written.
        proteome_protein_map (dict): A dictionary mapping protein_id to proteome labels.
        nolabel (str, optional): Label to use for proteins without a match. Default is "".
        sortlabels (bool, optional): Whether to sort the proteome labels. Default is False.
        uniq (bool, optional): Whether to skip identical lines in input. Default is False.

    Returns:
        int: The total number of clusters processed.
    """
    cluster_counter = -1
    prev_cluster = None
    prev_line = ""

    def _label_line(line):
        cluster_id, protein_id = line.rstrip("\n").split("\t")
        proteome_ids = proteome_protein_map.get(protein_id, nolabel)
        if sortlabels and "," in proteome_ids:
            proteome_ids = ",".join(sorted(set(proteome_ids.split(","))))
        return cluster_id, protein_id, proteome_ids

    with open(input_file, "r") as input_fh:
        with open(output_file, "w") as output_fh:
            output_fh.write(HEADER)
            if uniq:
                for line in input_fh:
                    if prev_line == line:
                        continue  # skip this line
                    else:
                        prev_line = line  # to check next
                    cluster_id, protein_id, proteome_ids = _label_line(line)
                    if cluster_id != prev_cluster:
                        prev_cluster = cluster_id
                        cluster_counter += 1
                        output_fh.write(
                            f"{cluster_counter}\t{protein_id}\t{proteome_ids}\t*\n"
                        )
                    else:
                        output_fh.write(
                            f"{cluster_counter}\t{protein_id}\t{proteome_ids}\t\n"
                        )
            else:
                for line in input_fh:
                    cluster_id, protein_id, proteome_ids = _label_line(line)
                    if cluster_id != prev_cluster:
                        prev_cluster = cluster_id
                        cluster_counter += 1
                        output_fh.write(
                            f"{cluster_counter}\t{protein_id}\t{proteome_ids}\t*\n"
                        )
                    else:
                        output_fh.write(
                            f"{cluster_counter}\t{protein_id}\t{proteome_ids}\t\n"
                        )

    return cluster_counter + 1


# functions for distributed approach
def re_label_proteins(
    input_file, output_file, proteome_protein_map, nolabel="", uniq=False
):
    """
    Relabel the input file using the provided mapping of protein IDs to proteome labels.
    If the input file already has labels, they will be preserved and updated with new labels.

    Args:
    - input_file (str): Path to the input file.
    - output_file (str): Path where the relabeled file will be written.
    - proteome_protein_map (dict): A mapping from protein IDs to proteome labels.
    - nolabel (str): Default label to use if no label is found for a protein.
    - uniq (bool): If True, duplicate lines in the input file are skipped.

    Returns:
    - None: This function modifies the output file in-place and does not return any value.

    Note:
        cluster_id (first column) will be converted to cluster_counter when combining all chunks.
    """
    already_labelled = False
    with open(input_file, "r") as input_fh:
        # check if this file has been labelled before
        first_line = input_fh.readline()
        fields_count = len(first_line.split("\t"))
        if fields_count == 3:  # previously labelled (3 columns file)
            already_labelled = True
        elif fields_count == 2:  # never labelled (2 columns file)
            pass
        else:
            raise ValueError(
                f"    => ERROR: file '{input_file}' has unexpected format: number of columns is {fields_count}: {first_line}"
            )

        input_fh.seek(0)  # reset to beginning of file to include first line
        with open(output_file, "w") as output_fh:
            if already_labelled:  # previously labelled
                for line in input_fh:
                    cluster_id, protein_id, old_proteome_ids = line.rstrip("\n").split(
                        "\t"
                    )
                    combined_proteome_ids = old_proteome_ids
                    new_proteome_ids = proteome_protein_map.get(protein_id, nolabel)
                    if combined_proteome_ids:
                        if new_proteome_ids:  # if there is data to combine
                            combined_proteome_ids += "," + new_proteome_ids  # append
                    else:
                        combined_proteome_ids = new_proteome_ids
                    output_fh.write(
                        f"{cluster_id}\t{protein_id}\t{combined_proteome_ids}\n"
                    )
            else:  # never labelled
                if uniq:
                    prev_line = ""
                    for line in input_fh:
                        if prev_line == line:
                            continue  # skip this line
                        else:
                            prev_line = line  # to check next
                        cluster_id, protein_id = line.rstrip("\n").split("\t")
                        output_fh.write(
                            f"{cluster_id}\t{protein_id}\t{proteome_protein_map.get(protein_id, nolabel)}\n"
                        )
                else:
                    for line in input_fh:
                        cluster_id, protein_id = line.rstrip("\n").split("\t")
                        output_fh.write(
                            f"{cluster_id}\t{protein_id}\t{proteome_protein_map.get(protein_id, nolabel)}\n"
                        )


def _find_next_newline(fp, filesize, startpos):
    """
    Return the position of the first newline found after the given starting position in the file.

    Args:
    - fp (file object): The file object to read from.
    - filesize (int): The size of the file in bytes.
    - startpos (int): The position in the file from where the search for the next newline starts.

    Returns:
    - int: The position (in bytes) of the first newline character after the starting position, or the file size if no newline is found.
    """
    if startpos > filesize:
        return filesize
    if startpos <= 0:
        return 0
    fp.seek(startpos - 1)
    fp.readline()  # to avoid truncatedlines
    return fp.tell()


def clone_file(filename, copies_num, outprefix, buffersize=BUFFERSIZE):
    """
    create multiple copies of an input file by writing each line to multiple output files simultaneously

    Args:
        filename (str): Path to the input file.
        copies_num (int): Number of copies to create.
        outprefix (str): Prefix for the output file names.
        buffersize (int): Number of bytes to read and write at a time. Default is 1MB.

    Returns:
        list: A list of filenames for the cloned files.
    """
    cloned_files = [outprefix + str(i) for i in range(copies_num)]
    file_handles = [open(filename, "wb") for filename in cloned_files]

    try:
        with open(filename, "rb") as infile:
            while True:
                chunk = infile.read(buffersize)
                if not chunk:
                    break  # end of file
                for handle in file_handles:
                    handle.write(chunk)
    finally:
        for handle in file_handles:
            handle.close()

    return cloned_files


def simply_label(filename, my_proteome_protein_map, nolabel="", uniq=False):
    """
    Annotates identifiers in a file using a provided mapping and writes the updated contents to a temporary file. The temporary file is then renamed to the input filename.

    Args:
        filename (str): The path to the file to be annotated.
        my_proteome_protein_map (dict): A dictionary mapping protein identifiers to labels.
        nolabel (str, optional): A string to mark when no mapping was found. Default is an empty string.
        uniq (bool, optional): A flag indicating whether to unique input lines. Default is False.

    Returns:
        None: The file is overwritten with the labeled data.
    """
    labelled_file = filename + "_tmp"
    re_label_proteins(
        filename, labelled_file, my_proteome_protein_map, nolabel=nolabel, uniq=uniq
    )
    os.rename(labelled_file, filename)


def lock_and_label(
    file_chunk,
    my_proteome_protein_map,
    nolabel="",
    uniq=False,
    timeout=90,
    workerid=None,
):
    """
    Annotates identifiers in a chunk of a file while ensuring that no other processes concurrently modify the file.
    Uses a file lock mechanism to prevent race conditions during the annotation process.

    Args:
        file_chunk (str): The path to the chunk of the file to be annotated.
        my_proteome_protein_map (dict): A dictionary mapping protein identifiers to labels.
        nolabel (str, optional): A string to mark when no mapping was found. Default is an empty string.
        uniq (bool, optional): A flag indicating whether to unique input lines. Default is False.

    Returns:
        int: Returns 1 if labeling was successful, 0 if it timed out.
    """
    lock = FileLock(file_chunk + ".lock")
    my_timeout = timeout + randint(0, 20)
    try:
        with lock.acquire(timeout=my_timeout):
            simply_label(
                file_chunk, my_proteome_protein_map, nolabel=nolabel, uniq=uniq
            )
            return 1
    except Timeout:
        if workerid is None:
            eprint(
                f"    ERROR: giving up waiting to acquire file lock for 'f{file_chunk}' after trying for {secs2time(my_timeout)}"
            )
        else:
            eprint(
                f" [{workerid}] ERROR: giving up waiting to acquire file lock for 'f{file_chunk}' after trying for {secs2time(my_timeout)}"
            )
        return 0


def combine_output_chunks(output_chunk_paths, out_file, sortlabels=False, uniq=False):
    """
    Combines all file chunks into a single output file, assigning unique sequential cluster identifiers.
    Optionally, it sorts and deduplicates labels in each line.

    Args:
        output_chunk_paths (list of str): List of file paths to the output chunks.
        out_file (str): Path to the combined output file.
        sortlabels (bool, optional): If True, sorts and deduplicates labels in each line.
        uniq (bool, optional): If True, removes consecutive duplicate lines. Default is False.

    Returns:
        int: The total number of unique clusters written to the output file.
    """
    cluster_counter = -1
    prev_cluster = None
    prev_line = ""
    with open(out_file, "w") as output_fh:
        output_fh.write(HEADER)
        for chunk_file in output_chunk_paths:
            with open(chunk_file, "r") as input_fh:
                for line in input_fh:
                    if uniq:
                        if prev_line == line:
                            continue  # skip this line
                        else:
                            prev_line = line  # to check next
                    cluster_id, protein_id, proteome_ids = line.rstrip("\n").split("\t")
                    if sortlabels and "," in proteome_ids:
                        proteome_ids = ",".join(sorted(set(proteome_ids.split(","))))
                    if cluster_id != prev_cluster:
                        prev_cluster = cluster_id
                        cluster_counter += 1
                        output_fh.write(
                            f"{cluster_counter}\t{protein_id}\t{proteome_ids}\t*\n"
                        )
                    else:
                        output_fh.write(
                            f"{cluster_counter}\t{protein_id}\t{proteome_ids}\t\n"
                        )

    return cluster_counter + 1


# for nochunks:
def combine_output_files(output_file_paths, out_file, sortlabels=False, uniq=False):
    """
    Combines multiple output files into a single file, assigning unique sequential cluster identifiers.
    Optionally, it sorts and deduplicates labels in each line.

    Args:
        output_file_paths (list): List of file paths to the files to be combined.
        out_file (str): Path to the final combined output file.
        sortlabels (bool, optional): If True, sorts and deduplicates the labels in each line.
                                     Default is False.
        uniq (bool, optional): If True, removes consecutive duplicate lines present in input.
                               Default is False.

    Returns:
        int: The total number of unique clusters written to the output file.
    """
    cluster_counter = -1
    prev_cluster = None
    prev_line = ""
    input_filehandles = [open(filename, "r") for filename in output_file_paths]

    try:
        with open(out_file, "w") as output_fh:
            output_fh.write(HEADER)
            for lines in zip(
                *input_filehandles
            ):  # simultaneously go through all file line by line
                cluster_id, protein_id, proteome_ids = (
                    lines[0].rstrip("\n").split("\t")
                )  # take all info from the first file
                if uniq:
                    if prev_line == "\t".join([cluster_id, protein_id]):
                        continue  # skip this line
                    else:
                        prev_line = "\t".join([cluster_id, protein_id])  # to check next

                for line in lines[1:]:
                    new_proteome_ids = line.rstrip("\n").split("\t")[
                        2
                    ]  # take only proteome_id from all other files
                    if new_proteome_ids:
                        proteome_ids += "," + new_proteome_ids

                proteome_ids = proteome_ids.lstrip(
                    ","
                )  # remove initial space if present
                if sortlabels and "," in proteome_ids:
                    proteome_ids = ",".join(sorted(set(proteome_ids.split(","))))

                if cluster_id != prev_cluster:
                    prev_cluster = cluster_id
                    cluster_counter += 1
                    output_fh.write(
                        f"{cluster_counter}\t{protein_id}\t{proteome_ids}\t*\n"
                    )
                else:
                    output_fh.write(
                        f"{cluster_counter}\t{protein_id}\t{proteome_ids}\t\n"
                    )

    finally:  # close opened file handles
        for file in input_filehandles:
            file.close()

    return cluster_counter + 1


def initializer(file_chunks_arg, args_arg):
    """
    initializer to set global variables for workers
    """
    global file_chunks, chunks_count, args
    file_chunks = file_chunks_arg
    chunks_count = len(file_chunks)
    args = args_arg


def worker_process(batch_files):
    """
    process to create a partial dictionary from a batch of files and annotate all file chunks
    """
    workerid = int(current_process().name.split("-")[1]) - 1  # 0..threads-1
    # eprint(f"   [{workerid}]: reading protein_ids from these files: {batch_files}") #debug

    # load partial dictionary for assigned batch of files
    my_proteome_protein_map = create_proteome_protein_map(
        batch_files, args, read_method=READMETHOD
    )
    # eprint(f"   [{workerid}]: created map with {len(my_proteome_protein_map)} entries") #debug

    if args.nochunks:
        my_file = file_chunks[workerid]  # each worker works on only one (big) file
        # eprint(f"  [{workerid}] processing my_file {my_file}")) #debug
        simply_label(
            my_file, my_proteome_protein_map, nolabel=args.nolabel, uniq=args.uniq
        )
    else:
        # dictionary to track labeling status for each file chunk
        chunk_status = {file_chunk: False for file_chunk in file_chunks}

        while not all(chunk_status.values()):  # repeat until all chunks are labelled
            still_to_label = [
                chunk for chunk, labeled in chunk_status.items() if not labeled
            ]
            for i in sample(range(len(still_to_label)), len(still_to_label)):
                file_chunk = still_to_label[i]
                # eprint(f"   [{workerid}]:   now labelling file chunk {file_chunk}") #debug
                result = lock_and_label(
                    file_chunk,
                    my_proteome_protein_map,
                    nolabel=args.nolabel,
                    uniq=args.uniq,
                )
                if result == 1:
                    chunk_status[file_chunk] = True
                    # eprint(f"  [{workerid}] acquired lock on {file_chunk} and labelled") #debug
                else:
                    eprint(
                        f"  [{workerid}] giving up on acquiring lock on {file_chunk}, will retry later. Still to work on: {[chunk for chunk, labeled in chunk_status.items() if not labeled]}"
                    )  # debug

        # process chunks randomly to avoid workers fighting for the same file_chunk
        # for i in sample(range(chunks_count), chunks_count):
        #    file_chunk = file_chunks[i]
        #    # eprint(f"   [{workerid}]:   now labelling file chunk #{i}") #debug
        #    lock_and_label(file_chunk, my_proteome_protein_map, nolabel=args.nolabel, uniq=args.uniq)
        # eprint(f"   [{workerid}]: labelled {chunks_count} files") #debug

    return len(my_proteome_protein_map)


if __name__ == "__main__":
    # ===============
    # 0: argument parsing and input/output files checking
    # ===============
    initial_secs = time.time()  # for total time count
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
    eprint(f" .-- BEGUN {timestamp} --.")

    args, fasta_files = check_args(DESCRIPTION)

    if args.distributed:
        # ===============
        # P1: setup: assign input file and fasta files among workers in batches
        # ===============
        start_secs = time.time()
        eprint(f" |-- ...")
        if args.nochunks:
            if args.threads == 1:
                # make a single copy of the input file if threads==1
                output_chunk_paths = [args.input_file + "_0"]
                shutil.copy(args.input_file, output_chunk_paths[0])
            else:
                # make args.threads cloned copies of the input file
                output_chunk_paths = clone_file(
                    args.input_file, args.threads, args.input_file + "_"
                )
                eprint(f" |-- input file cloned to each of the {args.threads} workers")
        else:
            output_chunk_paths, _ = split_file(
                args.input_file, args.chunksize, args.input_file + "_"
            )
            eprint(f" |-- input file split into {len(output_chunk_paths)} chunks")

        # clean up any temporary files leftover from previous executions, if any
        tmp_files = [filename + "_tmp" for filename in output_chunk_paths] + [
            filename + ".lock" for filename in output_chunk_paths
        ]
        # eprint(f" |-- deleting tmp files {tmp_files}") #debug
        delete_files(tmp_files)

        # split fasta_files into batches of batchsize files for each worker
        batches = [
            fasta_files[i : i + args.batchsize]
            for i in range(0, len(fasta_files), args.batchsize)
        ]
        eprint(
            f" |-- assigned fasta files in {len(batches)} batches of {args.batchsize} files to {min(args.threads, len(batches))} worker{'s' if min(args.threads, len(batches)) > 1 else ''}"
        )

        eprint(
            " |-- workers setup {} -- Elapsed: {} --".format(
                time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),
                elapsed_time(start_secs),
            )
        )

        # ===============
        # P2: run worker processes to label all files in a distributed fashion
        # ===============
        start_secs = time.time()
        processed_ids_counts = []
        with Pool(
            args.threads, initializer=initializer, initargs=(output_chunk_paths, args)
        ) as pool:
            if args.progress:
                for processed_ids_count in tqdm(
                    pool.imap_unordered(worker_process, batches), total=len(batches)
                ):
                    processed_ids_counts.append(processed_ids_count)

            else:
                for processed_ids_count in pool.imap_unordered(worker_process, batches):
                    processed_ids_counts.append(processed_ids_count)

        total_identifiers_processed = sum(processed_ids_counts)
        eprint(
            " |-- labelling completed {} -- Elapsed: {}, {} ids/s --".format(
                time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),
                *elapsed_time(start_secs, total_identifiers_processed),
            )
        )
        eprint(
            f" |-- processed {total_identifiers_processed} protein ids from {len(fasta_files)} proteome files"
        )

        # ===============
        # P3: combine the labelled files assigning cluster ids
        # ===============
        start_secs = time.time()
        if args.nochunks:
            # combine the complete files edited independently by the workers, assigning sequential cluster identifiers
            clusters_count = combine_output_files(
                output_chunk_paths,
                args.out_file,
                sortlabels=args.sortlabels,
                uniq=args.uniq,
            )
        else:
            # concatenate the chunks edited in parallel by the workers, assigning sequential cluster identifiers
            clusters_count = combine_output_chunks(
                output_chunk_paths,
                args.out_file,
                sortlabels=args.sortlabels,
                uniq=args.uniq,
            )
        eprint(f" |-- processed {clusters_count} clusters")
        eprint(
            " |-- final file created {} -- Elapsed: {}, {} clusters/s --".format(
                time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),
                *elapsed_time(start_secs, clusters_count),
            )
        )

        # clean up all temporary files
        tmp_files = (
            output_chunk_paths
            + [filename + "_tmp" for filename in output_chunk_paths]
            + [filename + ".lock" for filename in output_chunk_paths]
        )
        # eprint(f" |-- deleting tmp files {tmp_files}") #debug
        delete_files(tmp_files)
    else:  # sequential approach
        # ===============
        # S1: create protein->proteome map
        # ===============
        start_secs = time.time()
        eprint(f" |-- ...")
        proteome_protein_map = create_proteome_protein_map(
            fasta_files, args, read_method=READMETHOD
        )
        eprint(
            f" |-- processed {len(proteome_protein_map)} protein ids from {len(fasta_files)} proteome files"
        )
        eprint(
            " |-- mapping completed {} -- Elapsed: {}, {} ids/s --".format(
                time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),
                *elapsed_time(start_secs, len(proteome_protein_map)),
            )
        )

        # ===============
        # S2: assign proteome labels to protein identifiers
        # ===============
        start_secs = time.time()
        clusters_count = label_proteins(
            args.input_file,
            args.out_file,
            proteome_protein_map,
            nolabel=args.nolabel,
            sortlabels=args.sortlabels,
            uniq=args.uniq,
        )
        eprint(f" |-- processed {clusters_count} clusters")
        eprint(
            " |-- labelling completed -- Elapsed: {}, {} clusters/s --".format(
                *elapsed_time(start_secs, clusters_count),
            )
        )

    # total time
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
    eprint(
        " |-- total time: {} --".format(
            elapsed_time(initial_secs),
        )
    )
    eprint(f" '-- ENDED {timestamp} --'")
