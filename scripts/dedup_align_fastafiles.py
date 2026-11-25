#!/usr/bin/env python3
# changelog
# Sat 28 Dec 22:30:18 GMT 2024 0.1 started
# Sat 28 Dec 23:32:22 GMT 2024 0.9 coded
# Sun 29 Dec 09:01:28 GMT 2024 1.0 bugfixes
# Mon 20 Oct 11:32:09 BST 2025 1.1 added -m to create mapping files
# Tue 17 Feb 10:23:23 GMT 2026 1.2 added --all and --consensus

# imports
import os
from collections import defaultdict
from Bio import SeqIO, AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord

import hashlib
import argparse
from multiprocessing import Pool, current_process, cpu_count
import subprocess
import glob
import time
import sys
from tqdm import tqdm
from random import shuffle

# constants
DESCRIPTION = """
Wrapper around Clustal Omega that can find and deduplicate sequences across FASTA files, run alignments and optionally create mapping files of duplicates.

 * Accepts either a directory of FASTA files or a text file listing FASTA paths (mutually exclusive; one required)
 * Deduplicates identical sequences across input FASTA files (optionally recording mapping files of duplicated sequence occurrences (-m))
 * Alternatively align all sequences (-a), skipping deduplication step
 * Runs Clustal Omega to produce multiple sequence alignments for each input FASTA, using a configurable number of threads
 * Writes aligned output files into a specified output directory (created if missing)
 * Tweakable parallelism: a global worker pool (--threads) manages file-level parallelism; Clustal Omega itself receives a separate thread count (--align_threads) per job
 * Optionally write the Clustal consensus line in a separate file for each alignment

Example calls:
    dedup_align_fastafiles.py -d fasta_dir -o aligned_out -t 4 -at 2 # align all FASTAs in a directory with 4 workers and 2 threads per Clustal Omega run
    dedup_align_fastafiles.py -l my_fastas.txt -o aligned_out -m -q # align files listed in a text file, create duplicate mapping files, show progress
    dedup_align_fastafiles.py -l my_fastas.txt -a -o aligned_out # align all sequences of listed files, skipping deduplication step
"""
CLUSTALO_EXE = "clustalo"  # Path or name of Clustal Omega executable


# helper functions
def secs2time(secs):
    """Converts a time duration in seconds to a human-readable format (hours, minutes, seconds)."""
    minutes, seconds = divmod(secs, 60)
    hours, minutes = divmod(minutes, 60)
    return "{:02.0f}h {:02.0f}m {:02.0f}s".format(hours, minutes, seconds)


def elapsed_time(start_time, work_done=None, invert=False):
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
    if invert:
        process_speed = round(process_time / work_done, 2)
    else:
        process_speed = round(work_done / process_time, 2)
    return secs2time(process_time), process_speed


def exit_with_error(message: str, code: int = 1):
    """Prints an error message to stderr and exits the program with the specified exit code."""
    eprint(f" => {message}")
    sys.exit(code)


def eprint(*myargs, **kwargs):
    """Prints the provided arguments to stderr."""
    print(*myargs, file=sys.stderr, **kwargs)


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
            raise argparse.ArgumentTypeError("{} is not an integer".format(value))
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
            Print custom text before the default help message.
            """
            print(DESCRIPTION)
            super().print_help(*args, **kwargs)

    parser = CustomArgumentParser(
        description="Clustal Omega wrapper with sequence deduplication"
    )

    # General input options
    input_group = parser.add_argument_group("Input Options (choose one)")
    exclusive_group = input_group.add_mutually_exclusive_group(required=True)
    exclusive_group.add_argument(
        "-d",
        "--fasta_dir",
        type=str,
        required=False,
        dest="input_fasta_dir",
        help="Directory containing fasta files to be aligned.",
    )
    exclusive_group.add_argument(
        "-l",
        "--fasta_list",
        type=is_valid_file,
        required=False,
        dest="input_fasta_list",
        help="Path to a text file containing a list of FASTA files.",
    )
    exclusive_group.add_argument(
        "-s",
        "--single_file",
        type=str,
        nargs="+",
        required=False,
        dest="fasta_file",
        help="Fasta file name(s), space separated",
    )

    # General options
    parser.add_argument(
        "-t",
        "--threads",
        type=positive_integer,
        default=1,
        help="Number of threads for parallel processing.",
    )
    parser.add_argument(
        "-e",
        "--extension",
        type=str,
        default=".fa",
        help="File extension for FASTA files (default: .fa).",
    )

    parser.add_argument(
        "-a",
        "--all",
        action="store_true",
        required=False,
        help="Skip the deduplication step and align all sequences",
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
        "-m",
        "--mapfiles",
        action="store_true",
        required=False,
        help="Create mapping files with duplicated sequences",
        default=False,
    )
    parser.add_argument(
        "-c",
        "--consensus",
        action="store_true",
        required=False,
        help="Writes consensus line to file .cns, if present (clustal format only)",
        default=False,
    )

    clustalo_group = parser.add_argument_group("Clustal Omega specific arguments")

    clustalo_group.add_argument(
        "-o",
        "--out_dir",
        type=str,
        required=True,
        help="Output directory for aligned files.",
    )

    clustalo_group.add_argument(
        "-at",
        "--align_threads",
        type=positive_integer,
        default=1,
        help="Number of threads for Clustal Omega.",
    )

    clustalo_group.add_argument(
        "-f",
        "--force",
        action="store_true",
        help="Force overwrite of existing output files.",
    )

    clustalo_group.add_argument(
        "-st",
        "--seqtype",
        type=str,
        choices=["Protein", "DNA"],
        default="Protein",
        help="Specify the sequence type (e.g., 'Protein' or 'DNA').",
    )

    clustalo_group.add_argument(
        "-of",
        "--outfmt",
        type=str,
        choices=["fasta", "clustal", "stockholm"],
        default="fasta",
        help="Specify the MSA output format.",
    )

    args = parser.parse_args()

    # Validate input arguments
    if args.input_fasta_dir and not os.path.isdir(args.input_fasta_dir):
        exit_with_error(f"ERROR: No such directory '{args.input_fasta_dir}'", 2)

    if not os.path.isdir(args.out_dir):
        try:
            os.makedirs(args.out_dir)
            eprint(f" |-- Created output directory: {args.out_dir}")
        except PermissionError:
            exit_with_error(f"ERROR: Cannot create directory '{args.out_dir}'", 1)

    eprint(f" |-- Output directory: {args.out_dir}")

    if args.all and args.mapfiles:
        exit_with_error(
            f"ERROR: Cannot create mapping file when deduplication step is skipped", 22
        )

    if args.consensus and args.outfmt != "clustal":
        exit_with_error(
            f"ERROR: Cannot create consensus file when alignment format is {args.outfmt}",
            22,
        )

    if args.all:
        eprint(f" |-- Aligning all sequences, skipping deduplication")
    if args.mapfiles:
        eprint(f" |-- Creating mapping files")
    if args.consensus:
        eprint(f" |-- Creating consensus files")

    if args.threads > cpu_count():
        args.threads = cpu_count()
        eprint(f" |-- WARNING: only {args.threads} threads available")

    return args


def deduplicate_fasta(input_fasta, deduplicated_fasta, mapping_file=None):
    """
    Deduplicates sequences in a FASTA file, writing unique sequences to a new file
    and a mapping of which sequences were squashed together.
    """
    seq_to_headers = defaultdict(list)

    with open(deduplicated_fasta, "w") as dedup_fasta:
        for record in SeqIO.parse(input_fasta, "fasta"):
            sequence = str(record.seq).upper()
            sequence_hash = hashlib.sha1(sequence.encode()).hexdigest()  # Compute hash
            if sequence_hash not in seq_to_headers:  # unique sequence found
                dedup_fasta.write(f">{record.id}\n{sequence}\n")
            seq_to_headers[sequence_hash].append(record.id)

    # Optionally write a mapping file
    if mapping_file is not None:
        with open(mapping_file, "w") as map_file:
            for headers in seq_to_headers.values():
                map_file.write(f"{' '.join(headers)}\n")

    return seq_to_headers


def align_sequences(input_file, output_file, args):
    """
    Aligns sequences using Clustal Omega.
    """
    clustalo_call = [
        CLUSTALO_EXE,
        "--infile",
        input_file,
        "--outfile",
        output_file,
        "--outfmt",
        args.outfmt,
        "--threads",
        str(args.align_threads),
        "--seqtype",
        args.seqtype,
    ]
    if os.path.isfile(output_file) and args.force:
        clustalo_call.append("--force")

    # eprint(f"\nRunning Clustal Omega with command:\n {' '.join(clustalo_call)}") #debug
    try:
        subprocess.run(clustalo_call, check=True)
        return 1  # success run

    except subprocess.CalledProcessError as e:
        eprint(f"Error during alignment of {input_file}")
        eprint(f"Command: {' '.join(clustalo_call)}")
        eprint(f"Exit Code: {e.returncode}")
        eprint(f"Error Output: {e.stderr}")
        return 0  # failure


def expand_alignment(
    aligned_fasta,
    seq_to_headers,
    final_fasta,
    alignment_format="fasta",
    consensus_file=None,
):
    """
    Expands the alignment by duplicating sequences based on the mapping.
    Preserves clustal consensus if present, optionally writes it to file.

    NOTE: long sequence headers are cropped!
    """
    alignment = AlignIO.read(aligned_fasta, alignment_format)
    if (
        consensus_file is not None
        and "clustal_consensus" in alignment.column_annotations
    ):
        with open(consensus_file, "w") as consensusfh:
            consensusfh.write(
                f"'{alignment.column_annotations['clustal_consensus']}'\n"
            )

    record_map = {rec.id: rec for rec in alignment}
    expanded_records = []
    for seq_hash, headers in seq_to_headers.items():
        original = record_map[headers[0]]

        for header in headers:
            expanded_records.append(
                SeqRecord(seq=original.seq, id=header, description="")
            )
    expanded_alignment = MultipleSeqAlignment(expanded_records)
    expanded_alignment.column_annotations = alignment.column_annotations.copy()

    AlignIO.write(expanded_alignment, final_fasta, alignment_format)


def old_expand_alignment(
    aligned_fasta, seq_to_headers, final_fasta, alignment_format="fasta"
):
    """
    Expands the alignment by duplicating sequences based on the mapping.

    NOTE: clustal consensus line is lost
    """
    aligned_records = SeqIO.to_dict(SeqIO.parse(aligned_fasta, alignment_format))

    with open(final_fasta, "w") as final_fasta_file:
        for seq_hash, headers in seq_to_headers.items():
            aligned_seq = aligned_records[headers[0]].seq
            for header in headers:
                final_fasta_file.write(f">{header}\n{aligned_seq}\n")

    if alignment_format != "fasta":  # convert to desired format
        records = list(SeqIO.parse(final_fasta, "fasta"))
        SeqIO.write(records, final_fasta, alignment_format)


def initializer(args_arg):
    """
    initializer to set global variables for workers
    """
    global args
    args = args_arg


def worker_process(input_file):
    """process to run for multithreading operation"""
    workerid = int(current_process().name.split("-")[1]) - 1  # Worker ID for debugging
    # eprint(f"[Worker {workerid}/{args.threads}] Processing {my_file}") #debug

    if not os.path.isfile(input_file):
        eprint(f" |-- ERROR: no such file {input_file}")
        return 0

    # files involved
    base_name = os.path.splitext(os.path.basename(input_file))[0]
    output_file = os.path.join(args.out_dir, f"{base_name}.aln")
    deduplicated_fasta = os.path.join(args.out_dir, f"{base_name}.ded.tmp")
    aligned_fasta = os.path.join(args.out_dir, f"{base_name}.aln.tmp")
    if args.mapfiles:  # optional creation of mapping file
        mapping_file = os.path.join(args.out_dir, f"{base_name}.map")
    else:
        mapping_file = None
    if args.consensus:  # optional creation of consensus file
        consensus_file = os.path.join(args.out_dir, f"{base_name}.cns")
    else:
        consensus_file = None

    if os.path.isfile(output_file):
        if not args.force:
            eprint(
                f" |-- WARNING: Cowardly refusing to overwrite already existing file '{output_file}'. Use --force to force overwriting."
            )
            return 0

    if args.all:
        # align the sequences
        result = align_sequences(input_file, output_file, args)
        if result == 0:  # cleanup
            if os.path.isfile(output_file):
                os.remove(output_file)
        elif args.consensus:  # write consensus file
            alignment = AlignIO.read(output_file, "clustal")
            if "clustal_consensus" in alignment.column_annotations:
                with open(consensus_file, "w") as consensusfh:
                    consensusfh.write(
                        f"'{alignment.column_annotations['clustal_consensus']}'\n"
                    )
        return result
    else:
        # deduplicate the input FASTA
        seq_to_headers = deduplicate_fasta(input_file, deduplicated_fasta, mapping_file)

        if len(seq_to_headers) == 1:
            # won't align a single sequence
            eprint(
                f" |-- WARNING: '{input_file}' has a single unique sequence, skipping."
            )
            os.remove(deduplicated_fasta)
            return 0

        # align the deduplicated sequences
        result = align_sequences(deduplicated_fasta, aligned_fasta, args)

        if result == 0:
            # cleanup
            for filename in [deduplicated_fasta, aligned_fasta]:
                if os.path.isfile(filename):
                    os.remove(filename)
            return result

        # expand the alignment back to include original headers
        expand_alignment(
            aligned_fasta,
            seq_to_headers,
            output_file,
            alignment_format=args.outfmt,
            consensus_file=consensus_file,
        )

        # cleanup
        for filename in [deduplicated_fasta, aligned_fasta]:
            if os.path.isfile(filename):
                os.remove(filename)

    return result


if __name__ == "__main__":
    initial_secs = time.time()  # For total time count
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())

    eprint(f" .-- BEGUN {timestamp} --.")
    args = check_args(DESCRIPTION)

    # Determine input mode and validate
    if args.input_fasta_dir:
        eprint(f" |-- Input mode: Directory ({args.input_fasta_dir})")
        fasta_files = glob.glob(
            os.path.join(args.input_fasta_dir, f"*{args.extension}")
        )
        eprint(f" |-- Found {len(fasta_files)} FASTA files in directory.")
    elif args.input_fasta_list:
        eprint(f" |-- Input mode: File List ({args.input_fasta_list})")
        with open(args.input_fasta_list, "r") as f:
            fasta_files = [line.strip() for line in f if line.strip()]
        eprint(
            f" |-- Found {len(fasta_files)} FASTA files listed in {args.input_fasta_list}."
        )
    elif args.fasta_file:
        eprint(f" |-- Input mode: filenames ({args.fasta_file})")
        fasta_files = args.fasta_file
        for path in fasta_files:
            if not os.path.exists(path):
                exit_with_error(f"ERROR: No such file '{path}'", 2)
        eprint(f" |-- Operating on {len(fasta_files)} FASTA files")
    else:
        exit_with_error(
            "ERROR: No valid input specified. Please provide either --input_fasta_dir or --input_fasta_list or --single_file",
            1,
        )

    fasta_files_count = len(fasta_files)
    shuffle(
        fasta_files
    )  # if multiple scripts in parallel, each one will start working on different files
    results = []
    with Pool(
        processes=args.threads, initializer=initializer, initargs=(args,)
    ) as pool:
        if args.progress:
            for result in tqdm(
                pool.imap_unordered(worker_process, fasta_files),
                desc="aligning",
                total=fasta_files_count,
            ):
                results.append(result)
        else:
            for result in pool.imap_unordered(worker_process, fasta_files):
                results.append(result)

    successful_count = sum(results)

    # Only print stats if all files succeeded
    if successful_count != fasta_files_count:
        eprint(
            f" |-- ERROR: {fasta_files_count - successful_count} file(s) failed to process."
        )

    if successful_count == 0:
        eprint(" |-- ERROR: no alignment created")
    else:
        eprint(
            "|-- {} alignments completed -- Elapsed: {}, {} s/alignment --".format(
                successful_count,
                *elapsed_time(initial_secs + 1e-10, successful_count, invert=True),
            )
        )
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
    eprint(f" '-- ENDED {timestamp} --'")
