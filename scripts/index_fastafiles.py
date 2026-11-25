#!/usr/bin/env python3
# changelog
# Fri 29 Nov 2024 21:38:32 GMT 0.1 started
# Sat 30 Nov 2024 10:17:55 GMT 1.0 working singlethread
# Sat 30 Nov 2024 12:48:55 GMT 1.1 added multithreading

# imports
import os
import sys
import re
from glob import glob
import time
from sortedcontainers import SortedList
import argparse
from tqdm import tqdm
from typing import List, Optional
from ffdb import entry_generator, int_to_b64
from multiprocessing import Pool, current_process, cpu_count

# constants
FIELDSEPARATOR = "\t"  # delimiter for the index
BUFFERSIZE = 10 * 1048576  # 10Mb
DESCRIPTION = r"""
This script indexes FASTA files using ffdb (https://pypi.org/project/ffdb/), creates index files alongside the originals. It supports indexing a single file, a list of files or all the files with a given extension present in a directory.

 * Uses ffdb to build an index for each input FASTA file; output index filename is original_filename + ".idx".
 * Identifier extraction uses a regex pattern (default '^(.+|.*)$') to capture sequence IDs from FASTA headers
 * Processes files from one of three sources:
   - Single file: --input_file PATH
   - File-of-filenames: --list PATH (text file with one FASTA path per line)
   - Directory: --fasta_dir PATH with the possibility to filter by filename extension

Example calls:
    index_fastafiles.py --input_file mysequences.fasta
    index_fastafiles.py --list myfastafiles.tsv
    index_fastafiles.py --fasta_dir pigsequences/ -e .fa

The default pattern to capture fasta identifiers is '^(.+\|.*)$' but this can be modified with --pattern
E.g.
  --pattern '^(?:sp|tr)\|(\S*)\|\S* .*$' # to index by ACcession from UniProt fasta files
  --pattern '^(.+)\|.*$'               # to only get pid without other information after '|'

Parallel processing
    The files can be indexed in parallel by specifying --threads
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

    parser = CustomArgumentParser(description="Indexer for fasta files.")

    input_group = parser.add_mutually_exclusive_group(required=True)

    input_group.add_argument(
        "-i",
        "--input_file",
        dest="input_file",
        type=is_valid_file,
        help="Path to the fasta file to be indexed",
    )
    input_group.add_argument(
        "-l",
        "--list",
        dest="list",
        help="Path to a file containing a list of filepaths to index",
        type=is_valid_file,
    )

    input_group.add_argument(
        "-d",
        "--fasta_dir",
        dest="fasta_dir",
        type=str,
        help="Directory containing fasta files (with .fa extension unless -e specified)",
    )
    parser.add_argument(
        "-e",
        "--extension",
        type=str,
        required=False,
        default=".fa",
        help="Extension for files in fasta_dir; default '.fa'",
    )
    parser.add_argument(
        "-f",
        "--force",
        action="store_true",
        required=False,
        help="Force re-creation of existing index files",
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
    parser.add_argument(
        "-p",
        "--pattern",
        required=False,
        help=r"Regexp Capture pattern for identifiers to index in fasta header. Default is: '^(.+\|.*)$'",
        default=r"^(.+\|.*)$",
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=positive_integer,
        required=False,
        default=1,
        help="Number of threads for parallel processing",
    )

    args = parser.parse_args()

    if args.input_file is not None:
        eprint(f" |-- input_file: {args.input_file} will be indexed")
        fasta_files = [args.input_file]
    elif args.fasta_dir is not None:
        if not os.path.isdir(args.fasta_dir):
            exit_with_error(f"ERROR: No such directory '{args.fasta_dir}'", 2)
        fasta_files = glob(os.path.join(args.fasta_dir, "*" + args.extension))
        if not fasta_files:
            exit_with_error(
                f"ERROR: No matching '{args.extension}' files under '{args.fasta_dir}'.",
                2,
            )
        eprint(
            f" |-- fasta_dir: {args.fasta_dir} ({len(fasta_files)} files will be indexed)"
        )
    elif args.list is not None:
        fasta_files = []
        with open(args.list, "r") as fh:
            for line in fh.readlines():
                fasta_files.append(line.strip())
        eprint(
            f" |-- list of fasta filepaths: {args.list} ({len(fasta_files)} files will be indexed)"
        )
    else:
        exit_with_error(f"https://xkcd.com/2200/", 22)  # (we should never reach here)

    eprint(f" |-- identifier pattern: {args.pattern}")
    args.pattern = re.compile(args.pattern.encode("UTF-8"), re.MULTILINE)

    if args.threads > cpu_count():
        args.threads = cpu_count()
        eprint(f" |-- WARNING: only {args.threads} threads available")

    if args.threads > 1:
        eprint(f" |-- threads: {args.threads}")

    return args, fasta_files


def _format_indexes(ids, position, length):
    """
    return array of formatted indexes
    integer position and length of entry are converted to base64encoded number
    adapted from ffdb.py
    """
    indexes = []
    position_b64 = int_to_b64(position)
    length_b64 = int_to_b64(length)
    entryposition = f"{position_b64}-{length_b64}"
    for identifier in ids:
        indexes.append(f"{identifier}{FIELDSEPARATOR}{entryposition}\n")
    return indexes


def _find_patterns_in_entry(entry, pattern):
    """
    parse an entry, using regular expressions to capture desired fields
    adapted from ffdb.py's indexer.py
    """

    identifiers = list()
    match = pattern.search(entry)
    if match:
        for submatch in match.groups():
            # multiple capture patterns can be specified and we can index them all
            if submatch:
                # eprint(f"found {submatch}") #debug
                identifiers.append(submatch.decode())
    return identifiers


def _parse_ff_byentry(inputfile, args):
    """
    parse input file entry by entry
    adapted from ffdb.py's indexer.py
    """

    entries_count, indexes_count, skipped_count = (0, 0, 0)
    indexes = SortedList()
    entry = ""
    entrylength = 0
    outindex_filename = inputfile + ".idx"
    if not args.force and os.path.isfile(outindex_filename):
        eprint(
            f"Cowardly refusing to overwrite already existing index '{outindex_filename}'. Use --force to force overwriting."
        )  # debug
        return 0, 0, 0

    # init the entry generator
    gen = entry_generator(inputfile, ">", BUFFERSIZE)  # terminator: '>'
    # skip the first entry (consisting of a single '>')
    next(gen, None)  # None is the default return if generator is exhausted

    with open(outindex_filename, "w") as outindexfh:
        entryposition = (
            1  # position of first entry; since we removed the first '>' empty "entry"
        )
        for entry, entrylength in gen:
            # eprint(entry, entrylength) #debug
            entries_count += 1
            ids = _find_patterns_in_entry(entry, args.pattern)
            if ids:  # if identifiers found
                new_indexes = _format_indexes(ids, entryposition, entrylength)
                # eprint(f"collected new indexes for entry of size {entrylength}: {new_indexes}")
                indexes.update(new_indexes)
            else:  # we skip the entry since we found no identifiers
                skipped_count += 1
            entryposition += entrylength  # for next entry

        # now deal with last entry (since > is not a real terminator but an initiator)
        entries_count += 1
        with open(inputfile, "rb") as inputfh:
            inputfh.seek(entryposition)
            entry = inputfh.read()
            entrylength = len(entry)
            # eprint(f"read last entry: {entry}") #debug
            ids = _find_patterns_in_entry(entry, args.pattern)
            if ids:  # if identifiers found
                new_indexes = _format_indexes(
                    ids, entryposition, entrylength + 1
                )  # +1 to match the index produced by indexer.py (where we add a trailing '>' to input)
                # eprint(f"collected new indexes for last entry of size {entrylength}: {new_indexes}") #debug
                indexes.update(new_indexes)
            else:  # we skip the entry since we found no identifiers
                skipped_count += 1

        indexes_count = len(indexes)
        for index in indexes:
            outindexfh.write(index)

    return entries_count, indexes_count, skipped_count


def initializer(args_arg):
    """
    initializer to set global variables for workers
    """
    global args
    args = args_arg


def index_fasta_file(fasta_file):
    r"""
    uses routines adapted from ffdb.py to index a fasta file containing several sequences
    producing equivalent results to the following call of ffdb.py's indexer:
    indexer.py -e '>' -i '^(.+\|.*)$' -f <(cat $filepath; echo '>') > ${filepath}.idx
    """
    # workerid = int(current_process().name.split("-")[1]) - 1  # 0..threads-1
    # eprint(f" [{workerid}] working on {fasta_file}") #debug
    return _parse_ff_byentry(fasta_file, args)


if __name__ == "__main__":
    initial_secs = time.time()  # for total time count
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())

    eprint(f" .-- BEGUN {timestamp} --.")
    args, fasta_files = check_args()
    eprint(f" |...")

    results = []
    if args.progress:
        with Pool(
            processes=args.threads, initializer=initializer, initargs=(args,)
        ) as pool:
            for result in tqdm(
                pool.imap_unordered(index_fasta_file, fasta_files),
                total=len(fasta_files),
            ):
                results.append(
                    result
                )  # to ensure the results are consumed as they are produced

    else:
        with Pool(
            processes=args.threads, initializer=initializer, initargs=(args,)
        ) as pool:
            for result in pool.imap_unordered(index_fasta_file, fasta_files):
                results.append(result)
        # each worker returns entries_count, indexes_count, skipped_count

    # sum together the results coming from each worker
    transposed_results = zip(*results)
    total_entries, total_indexes, total_skipped = [
        sum(group) for group in transposed_results
    ]

    # stats
    if total_indexes > 0:
        eprint(
            f" |-- indexed {total_indexes} identifiers from {total_entries} fasta sequences"
        )
    else:
        eprint(f" |-- no sequence has been indexed")

    if total_skipped > 0:
        eprint(f" |-- {total_skipped} entries skipped")

    eprint(
        " |-- Elapsed: {}, parsed {} sequences/s --".format(
            *elapsed_time(initial_secs, total_entries),
        )
    )
    # total time
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
    eprint(f" '-- ENDED {timestamp} --'")
