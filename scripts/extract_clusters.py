#!/usr/bin/env python3
# changelog
# Tue 19 Nov 2024 19:16:32 GMT 1.0 coded
# Thu 21 Nov 2024 12:14:58 GMT 1.1 added --thread (multithreading) and --force (overwrite output)
# Thu 21 Nov 2024 20:06:51 GMT 1.2 added --all and its possible combination with --uniq
# Thu 12 Dec 2024 19:06:19 GMT 1.3 added --single_file
# Sun 15 Feb 23:35:44 GMT 2026 1.4 added --cluster
# Sun 16 Feb 22:04:17 GMT 2026 1.5 added --web for web retrieval, including sequence cache
# Wed 18 Feb 12:33:09 GMT 2026 1.6 added --uniqterse
# Wed 18 Feb 17:51:55 GMT 2026 1.7 count of duplicated sequences placed in header when using --uniqterse

# imports
import os
import re
import sys
import time
import argparse
from tqdm import tqdm
from hashlib import sha1
from typing import List, Optional
from ffdb import (
    split_file,
    calculate_blocksize,
    REES,
    b64_to_int,
    read_from_size,
    get_position_first,
)
from multiprocessing import Pool, current_process, cpu_count
from requests import get
from shutil import rmtree

# constants
HEADER = "cluster_id\tprotein_ids\tproteins_count\tproteomes_count\n"
REFASTAHEADER_STRIP = re.compile(r"^>.*?\n")  # fasta header removal
RESEQLEN_CLEAN = re.compile(r"\|[0-9]+$")  # clean seqlen info after identifier
BUFFERSIZE = 1048576  # 1Mb
DESCRIPTION = """
This script extracts FASTA sequences for protein identifiers grouped by cluster from a TSV (like the output of filter_clusters.py), retrieving sequences from ffdb-indexed FASTA files and writing one FASTA file per cluster.

 * Expects a TSV with at least two columns: cluster_id and protein_ids
   - Each entry in the protein_ids column is a space-separated token formatted as: proteome_id:protein_id
 * Writes, for each cluster, a FASTA file named <cluster_id>.fa containing the sequences of the proteins in that cluster
   - Sequences are extracted from FASTA files located in FASTA_DIR using ffdb index files (filename.idx)
   - Proteome IDs map to FASTA filenames; the script supports recovering filenames using optional prefixes/extensions


Sample input file:
    cluster_id\tprotein_ids\tproteins_count\tproteomes_count
    0\t35497:ENSSSCP00000055324|661 4698922:ENSSSCP00055011301|568 4698918:ENSSSCP00035021320|596\t3\t3

The script would create a fasta file called 0.fa containing the sequences listed for that cluster_id,
extracting them from the files corresponding to the proteome_id labels placed before each protein identifier.

Sample output file (sequence length, if present, is removed):
    >35497:ENSSSCP00000055324
    MSHESSQDRSSCRGSVVTNPNSIHEEDSVV[...]
    >4698922:ENSSSCP00055011301
    MNIKSLMKKSLVTCISFFFFFFSRNLVVRR[...]
    >4698918:ENSSSCP00035021320
    MNIKSLMKKSLVTCISFFFFFFSRNLVVRR[...]

Example calls:
    extract_clusters.py --input_file results_pig/Protein_clusters_m13.tsv --fasta_dir pig --out_dir clusters_pig -e .fa -p proteome_ -q #extract from a directory containing individual fasta files

    extract_clusters.py --input_file results_pig/Protein_clusters_m13.tsv --single_file pig --out_dir clusters_pig -e .fa -p proteome_ -q #extract from a single combined fasta file

Parallel processing (using --threads):
    The input file will be split in a number of chunks equal to the number of threads and each thread will work in parallel to create output cluster fasta files.

Options for comprehensive outputs:
    --uniq: to keep only unique sequences, merging identifiers into one header; e.g.:
    >4698918:ENSSSCP00035021316 4698918:ENSSSCP00035021317 4698918:ENSSSCP00035021319 4698918:ENSSSCP00035021320
    (protein identifiers sharing the same sequence are printed in the header, space separated)

    --all: to keep all the proteome tags from the input file into the output headers; e.g.:
    >35497,4229143,4698268,4698269,4698918,4698920,4698921,4698922,4698923,4698925,4698926:ENSSSCP00000055569
    (identifiers for proteomes which contain the protein are printed in the header, comma separated)

    The two options can be combined.

    NOTES:
        Pay attention that headers could become too long. These options may not be recommended
        if the number of proteins with the same sequence is expected to be very large.

        --all assumes that the input was produced using --all option of the filter_proteomes script
"""


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


# functions
def check_args():
    """
    parse arguments and check for error conditions
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

    parser = CustomArgumentParser(description="Clusters sequence extraction.")

    # Sequence origin options
    sequence_group = parser.add_argument_group("Fasta sequences from (choose one)")
    exclusive_group = sequence_group.add_mutually_exclusive_group(required=True)
    exclusive_group.add_argument(
        "-d",
        "--fasta_dir",
        type=str,
        required=False,
        dest="fasta_dir",
        help="Directory containing all proteome fasta files",
    )
    exclusive_group.add_argument(
        "-s",
        "--single_file",
        type=is_valid_file,
        required=False,
        dest="single_file",
        help="Path to a single fasta file where to extract all identifiers from",
    )
    exclusive_group.add_argument(
        "-w",
        "--web",
        type=str,
        required=False,
        dest="url",
        help="URL to web API where to fetch fasta sequences by UPI identifier e.g. 'https://rest.uniprot.org/uniparc/{UPI}.fasta' ({UPI} is where the UniParc identifier will be placed in the request",
    )

    # Input output options
    io_group = parser.add_argument_group("Input and Output")
    io_group.add_argument(
        "-i",
        "--input_file",
        type=is_valid_file,
        required=True,
        help="Path to the input file",
    )
    io_group.add_argument(
        "-o",
        "--out_dir",
        type=str,
        required=True,
        help="Directory where the cluster .fa files will be written",
    )

    # Comprehensive output options
    comprehensive_output_group = parser.add_argument_group("Comprehensive output")
    comprehensive_output_group.add_argument(
        "-u",
        "--uniq",
        action="store_true",
        required=False,
        help="Only keep unique sequences (merging the identifiers into one header) in case multiple protein identifiers have the same sequence.",
        default=False,
    )
    comprehensive_output_group.add_argument(
        "-a",
        "--all",
        action="store_true",
        dest="allproteomes",
        required=False,
        help="Keep all proteome labels in case the same protein_id is tagged to multiple proteomes in the input file.",
        default=False,
    )
    comprehensive_output_group.add_argument(
        "-ut",
        "--uniqterse",
        action="store_true",
        required=False,
        help="Only keep unique sequences (under the first seen identifier and without proteome label: shortest header)",
        default=False,
    )

    # Filtering options
    filter_group = parser.add_argument_group("Filtering options")
    filter_group.add_argument(
        "-m",
        "--minproteins",
        type=positive_integer,
        required=False,
        default=1,
        help="Optionally only extract clusters that contain at least this number of proteins. Note that this is checked first, before applying restricted filter.",
    )
    filter_group.add_argument(
        "-r",
        "--restrict",
        type=is_valid_file,
        required=False,
        help="Path to a file containing proteome identifiers. If given, only sequences belonging to the proteomes from that file will be extracted.",
    )
    filter_group.add_argument(
        "-c",
        "--cluster",
        type=str,
        required=False,
        nargs="+",
        help="Optionally specify which cluster(s) to extract, by cluster_id (multiple cluster_id can be specified, space separated)",
    )

    # Optional arguments:
    parser.add_argument(
        "-p",
        "--prefix",
        type=str,
        required=False,
        default="",
        help="Prefix for fasta filenames, to be added to proteome identifiers; e.g. 'proteome_'",
    )
    parser.add_argument(
        "-e",
        "--extension",
        type=str,
        required=False,
        default="",
        help="Extension for files in fasta_dir, to be added to proteome identifiers; e.g. '.fa'",
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=positive_integer,
        required=False,
        default=1,
        help="Number of threads for parallel processing.",
    )
    parser.add_argument(
        "-f",
        "--force",
        action="store_true",
        required=False,
        help="Force overwrite of existing cluster output files.",
        default=False,
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

    if args.single_file is not None:
        eprint(f" |-- single_file (for fasta sequences): {args.single_file}")
    elif args.fasta_dir:
        if not os.path.isdir(args.fasta_dir):
            exit_with_error(f"ERROR: No such directory '{args.fasta_dir}'", 2)
        eprint(f" |-- fasta_dir (for fasta sequences): {args.fasta_dir}")
    elif args.url:
        eprint(f" |-- web api url (for fasta sequences): {args.url}")
        args.sequence_cache = os.path.join(args.out_dir, "sequence_cache")
    else:  # should not happen
        exit_with_error(
            "ERROR: We've reached an unreachable state. Anything is possible. The limits were in our heads all along. Follow your dreams (xkcd#2200)",
            22,
        )

    if os.path.isdir(args.out_dir):
        try:
            temp_file = os.path.join(args.out_dir, "__TMPTESTFILE__")
            filehandle = open(temp_file, "w")
            delete_files([temp_file])
            if args.url is not None and not os.path.isdir(args.sequence_cache):
                os.mkdir(args.sequence_cache)
        except PermissionError:
            exit_with_error(f"ERROR: Cannot write into directory '{args.out_dir}'", 1)
    else:
        try:
            os.mkdir(args.out_dir)
            if args.url is not None and not os.path.isdir(args.sequence_cache):
                os.mkdir(args.sequence_cache)
        except PermissionError:
            exit_with_error(f"ERROR: Cannot create directory '{args.out_dir}'", 1)

    eprint(f" |-- out_dir: {args.out_dir}/")
    if args.cluster is not None:
        if len(args.cluster) == 1 and " " in args.cluster[0]:
            args.cluster = args.cluster[0].split(" ")
        eprint(f" |-- cluster(s) to extract: {' '.join(args.cluster)}")
    eprint(f" |-- min proteins threshold: {args.minproteins}")

    if args.threads > cpu_count():
        args.threads = cpu_count()
        eprint(f" |-- WARNING: only {args.threads} threads available")

    if args.threads > 1:
        args.chunksize = calculate_blocksize(args.input_file, args.threads)
        eprint(
            f" |-- threads: {args.threads}; input file will be split in chunks of max {args.chunksize} bytes"
        )

    if (args.uniq and args.uniqterse) or (args.allproteomes and args.uniqterse):
        exit_with_error(
            f"ERROR: option --uniqterse not compatible with --all or --uniq", 22
        )

    if args.uniq:
        eprint(f" |-- extracting only unique sequences, merging identifiers")
    elif args.uniqterse:
        eprint(
            f" |-- extracting only unique sequences, under the first seen identifier"
        )
    if args.allproteomes:
        eprint(f" |-- keeping all proteome tags in header")

    return args


def _extract_entry(identifier, filename, index_filename):
    """
    Extracts an entry from a flat-file database using ffdb.py.

    Args:
        entryid (str): The unique identifier of the entry to extract.
        filename (str): The filename of the database file.
        index_filename (str): The filename of the index file

    Returns:
        str or None: The extracted entry as a string, or None if not found.

    Raises:
        FileNotFoundError: If either the database or index file is not found.

    Reference:
        https://github.com/g-insana/ffdb.py

    """

    if not os.path.exists(filename):
        exit_with_error(f"ERROR: File '{filename}' does not exist.", 2)

    if not os.path.exists(index_filename):
        exit_with_error(f"ERROR: Index '{index_filename}' does not exist.", 2)

    # 1) extract position and entry size
    with open(index_filename, "r") as indexfh:
        position = get_position_first(indexfh, identifier)

    if position is None:
        eprint(f"    => WARNING: '{identifier}' not found in index; skipping")  # debug
        return None

    posmatch = REES.match(position)
    position, entry_length = posmatch.groups()

    # 2) decode position
    position = b64_to_int(position)
    entry_length = b64_to_int(entry_length)
    # eprint(f"  extracting '{identifier}' at position {position}, size {entry_length} from {filename}") #debug

    # 3) retrieve entry
    with open(filename, "r") as flatfilefh:
        entry = read_from_size(flatfilefh, position, entry_length)
    return entry


def retrieve_fasta_sequence(protein_id, weburl, sequence_cache):
    """
    Retrieves a fasta entry from an online api address

    Args:
        protein_id (str): The unique identifier of the fasta entry to extract
        weburl (str): URL to online service to retrieve sequnce by UPI
        sequence_cache (str): Path where to read/write sequence files

    Returns:
        sequence (str): without header or newline

    Note:
        it assumes {UPI} present in the weburl and UPI is present in protein_id
        as third field in protein_id format 'source_id|seqlen|UPI'
    """
    sequence = ""
    headers = {"user-agent": "extract_clusters/1.5", "Accept": "text/x-fasta"}
    parts = protein_id.split("|")
    if len(parts) >= 3:
        upi = parts[2]
        sequence = read_sequence_cache(upi, sequence_cache)
        if sequence == "":  # not cached
            url = weburl.replace("{UPI}", upi)
            # eprint(f"Fetching {upi} from weburl {url}") #debug
            r = get(url, headers=headers)
            if r.status_code in (200, 206):
                # eprint(f"got sequence {r.text}") #debug
                sequence = REFASTAHEADER_STRIP.sub("", r.text).rstrip()
                if sequence != "":
                    write_sequence_cache(upi, sequence, sequence_cache)
            else:
                eprint(
                    f"ERROR: problems retrieving {upi}, http status code is '{r.status_code}' ({r.reason})"
                )

    return sequence


def read_sequence_cache(upi, sequence_cache):
    """
    Retrieves a fasta sequence from cached file

    Args:
        upi (str): UniParc identifier of the sequence
        sequence_cache (str): Path whence to read the sequence file

    Returns:
        sequence (str): without header or newline
    """
    sequence = ""
    sequence_filename = os.path.join(sequence_cache, f"{upi}.fa")
    if os.path.isfile(sequence_filename):
        # eprint(f"Reading sequence for {upi} from cache") #debug
        with open(sequence_filename) as fh:
            sequence = fh.read()
    return sequence


def write_sequence_cache(upi, sequence, sequence_cache):
    """
    Writes a fasta sequence into a file in the cache

    Args:
        upi (str): UniParc identifier of the sequence
        squence (str): the sequence, without header
        sequence_cache (str): Path whence to read the sequence file

    Returns:
        Nothing
    """
    sequence_filename = os.path.join(sequence_cache, f"{upi}.fa")
    # eprint(f"Writing sequence for {upi} to cache") #debug
    with open(sequence_filename, "w") as fh:
        fh.write(sequence)


def extract_fasta_sequence(protein_id, filename):
    """
    Extracts a fasta entry from an indexed fasta file via ffdb.py

    Args:
        protein_id (str): The unique identifier of the fasta entry to extract.
        filename (str): The filename of the file from which to extract the entry.

    Returns:
        sequence (str): without header or newline

    Note:
        it assumes ffdb index to be located at filename + '.idx'

    """
    entry = _extract_entry(protein_id, filename, filename + ".idx")
    sequence = ""
    if entry is not None:
        sequence = REFASTAHEADER_STRIP.sub("", ">" + entry).rstrip()  # remove header
    return sequence


def extract_clusters(args):
    """
    Extracts sequences of proteins listed in input_file
    and creates new fasta files with them, one per cluster, under out_dir.

    Optionally, only those clusters which have at least a certain number of proteins in it will be
    extracted.
    Optionally, extract only proteins belonging to a given "restrict_set" of proteomes.
    Optionally, only extract specified cluster(s)
    Optionally, only print unique sequences, merging together duplicates.

    Args:
        args: argparse object

    Returns:
        int, int, int, int: total number of seen proteins, total number of seen clusters, total number of sequences and the total number of cluster files written.
    """
    input_file = args.input_file
    fasta_dir = args.fasta_dir
    weburl = args.url
    out_dir = args.out_dir
    uniq = args.uniq
    uniqterse = args.uniqterse
    restrict_set = args.restrict_set
    min_threshold = args.minproteins
    prefix = args.prefix
    extension = args.extension
    progress = args.progress
    allproteomes = args.allproteomes
    clusterslist = args.cluster

    protein_count = 0  # number of total proteins
    cluster_count = 0  # number of total clusters
    printed_sequences_count = (
        0  # number of total protein sequences printed in the output cluster files
    )
    created_clusterfiles_count = 0  # number of total output cluster files created

    input_file_size = os.path.getsize(input_file)
    if progress:
        pbar = tqdm(total=input_file_size, unit="B", position=args.pbar_position)

    with open(input_file, "r") as input_fh:
        for line in input_fh:  # for each cluster
            if progress:
                pbar.update(len(line))
            columns = line.rstrip("\n").split("\t", maxsplit=2)
            cluster_id, proteins_string = columns[0:2]
            if cluster_id == "cluster_id":
                # skip header
                continue
            cluster_count += 1
            protein2proteomes = {}
            for pair in proteins_string.split(" "):
                proteome_id, protein_id = pair.split(":", 1)
                if protein_id in protein2proteomes:
                    protein2proteomes[protein_id].append(proteome_id)
                else:
                    protein2proteomes[protein_id] = [proteome_id]

            # eprint(f"collected the following protein2proteomes matches for cluster #{cluster_id}: {protein2proteomes}") #debug
            this_line_proteins = len(protein2proteomes)
            protein_count += this_line_proteins
            if clusterslist is not None and cluster_id not in clusterslist:
                # not in the requested list
                # eprint(f"skipping cluster {cluster_id} as not requested") # debug
                continue
            if (
                this_line_proteins < min_threshold
            ):  # we skip if lower than the threshold
                # eprint(f"skipping cluster {cluster_id} as it has contains less proteins than {min_threshold}")  # debug
                continue

            # check if cluster file already present
            out_clusterfile = os.path.join(out_dir, cluster_id + ".fa")
            if not args.force and os.path.isfile(out_clusterfile):
                eprint(
                    f"Cowardly refusing to overwrite already existing file '{out_clusterfile}'. Use --force to force overwriting."
                )  # debug
                continue

            # reset dicts
            protein2sequence = {}
            sequencehash2protein = {}
            proteinid_synonyms = {}
            found_protein2proteome = {}

            # collect all sequences
            for protein_id, proteomes in protein2proteomes.items():
                # eprint(f"{protein_id} belongs to {proteomes}")
                proteome_id = proteomes[0]  # we will extract from the first one
                proteomes = set(proteomes)  # uniq in case of duplications
                if restrict_set is not None:
                    if not proteomes.intersection(restrict_set):
                        # eprint(f"skipping {protein_id} because none of the proteomes it belongs to are in the restrict_set") #debug
                        continue
                # eprint(f"extracting {protein_id} from {fasta_dir}/{prefix}{proteome_id}{extension}") #debug
                if args.single_file is not None:
                    sequence = extract_fasta_sequence(protein_id, args.single_file)
                elif fasta_dir is not None:  # fasta_dir
                    sequence = extract_fasta_sequence(
                        protein_id,
                        os.path.join(fasta_dir, prefix + proteome_id + extension),
                    )
                else:  # weburl, TODO
                    sequence = retrieve_fasta_sequence(
                        protein_id, weburl, args.sequence_cache
                    )
                if sequence == "":
                    eprint(
                        f"  WARNING skipping {proteome_id}:{protein_id} for cluster {cluster_id} because its sequence could not be retrieved"
                    )
                    continue

                protein_id = RESEQLEN_CLEAN.sub(
                    "", protein_id
                )  # clean seqlen info after identifier
                if allproteomes:
                    found_protein2proteome[protein_id] = ",".join(sorted(proteomes))
                else:  # only the first one
                    found_protein2proteome[protein_id] = proteome_id
                # keep only sequence from entry
                if (
                    uniq or uniqterse
                ):  # merge together those identifiers with same sequence
                    sequencehash = sha1(
                        sequence.replace("\n", "").upper().encode()
                    ).hexdigest()
                    same_seq_proteinid = sequencehash2protein.get(sequencehash, "")
                    if (
                        same_seq_proteinid != ""
                    ):  # if this sequence has been seen before
                        proteinid_synonyms[same_seq_proteinid] += (
                            " " + protein_id
                        )  # store a synonym
                    else:
                        sequencehash2protein[sequencehash] = (
                            protein_id  # store new sequence hash
                        )
                        protein2sequence[protein_id] = sequence
                        proteinid_synonyms[protein_id] = protein_id  # init
                else:
                    protein2sequence[protein_id] = sequence

            # print found sequences
            if len(protein2sequence):
                with open(out_clusterfile, "w") as output_fh:
                    created_clusterfiles_count += 1
                    for protein_id, sequence in protein2sequence.items():
                        printed_sequences_count += 1
                        if (
                            uniq
                        ):  # merged identifiers in header if multiple proteins with same sequence
                            identifiers = " ".join(
                                [
                                    found_protein2proteome[protein_id]
                                    + ":"
                                    + protein_id
                                    for protein_id in proteinid_synonyms.get(
                                        protein_id, protein_id
                                    ).split(" ")
                                ]
                            )
                            output_fh.write(f">{identifiers}\n{sequence}\n")
                        elif uniqterse:
                            if (
                                protein_id in proteinid_synonyms
                                and " " in proteinid_synonyms[protein_id]
                            ):
                                multiple = f"({len(proteinid_synonyms[protein_id].split(' '))})"
                            else:
                                multiple = ""
                            output_fh.write(f">{protein_id}{multiple}\n{sequence}\n")
                        else:
                            output_fh.write(
                                f">{found_protein2proteome[protein_id]}:{protein_id}\n{sequence}\n"
                            )
                # eprint(f"created {out_clusterfile} fasta file") #debug

    if progress:
        pbar.close()

    return (
        protein_count,
        cluster_count,
        printed_sequences_count,
        created_clusterfiles_count,
    )


def initializer(args_arg):
    """
    initializer to set global variables for workers
    """
    global args
    args = args_arg


def worker_process(my_file):
    """
    process to work in parallel extract_clusters on chunks of the input_file
    """
    workerid = int(current_process().name.split("-")[1]) - 1  # 0..threads-1
    args.input_file = my_file
    if args.progress:
        args.pbar_position = workerid
    # eprint(f" [{workerid}] working on {args.input_file}") #debug
    return extract_clusters(args)


if __name__ == "__main__":
    initial_secs = time.time()  # for total time count
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())

    eprint(f" .-- BEGUN {timestamp} --.")
    args = check_args()
    eprint(f" |...")

    args.restrict_set = None
    if args.restrict is not None:
        # read the list of proteomes we want to restrict output clusters to
        args.restrict_set = set()
        with open(args.restrict, "r") as input_fh:
            for line in input_fh:
                args.restrict_set.add(line.strip())

    if args.threads > 1:
        output_chunk_paths, _ = split_file(
            args.input_file, args.chunksize, args.input_file + "_"
        )
        with Pool(args.threads, initializer=initializer, initargs=(args,)) as pool:
            results = pool.map(worker_process, output_chunk_paths)

        # sum together the results coming from each worker
        transposed_results = zip(*results)
        total_proteins, total_clusters, printed_proteins, printed_clusters = [
            sum(group) for group in transposed_results
        ]

        # clean up of chunk files
        delete_files(output_chunk_paths)
    else:
        if args.progress:
            args.pbar_position = 0
        total_proteins, total_clusters, printed_proteins, printed_clusters = (
            extract_clusters(args)
        )

    eprint(
        f" |-- seen {total_proteins} protein_ids in {total_clusters} clusters (mean of {total_proteins/total_clusters:.2f} proteins/cluster)"
    )
    if printed_clusters:
        eprint(
            f" |-- extracted sequences for {printed_proteins} protein_ids in {printed_clusters} clusters (mean of {printed_proteins/printed_clusters:.2f} proteins/cluster)"
        )
        eprint(
            " |-- Elapsed: {}, parsed {} clusters/s --".format(
                *elapsed_time(initial_secs, total_clusters),
            )
        )
    else:
        eprint(f" |-- no cluster output files created")

    # cleanup sequence cache
    if args.url is not None and os.path.isdir(args.sequence_cache):
        # eprint(f"removing sequence cache {args.sequence_cache}") #debug
        rmtree(args.sequence_cache)
    # total time
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
    eprint(f" '-- ENDED {timestamp} --'")
