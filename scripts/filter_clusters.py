#!/usr/bin/env python3
# changelog
# Sat 16 Nov 2024 22:10:37 GMT 1.0 coded
# Sun 17 Nov 2024 17:52:43 GMT 1.1 added --progress
# Tue 19 Nov 2024 20:39:07 GMT 1.2 unique proteome_id:protein_id combinations only in output
# Sun  1 Dec 2024 14:07:05 GMT 1.3 added representative and seqlen_range columns
# Thu 19 Jun 11:26:28 BST 2025 1.4 added strict mode: ignore unlabelled protein_id (without # proteome_id)
# Thu 13 Nov 15:14:54 GMT 2025 1.5 allowed for more fields after seqlen
# Fri 12 Dec 21:58:17 GMT 2025 1.6 compute and print mode column (the most frequent protein length among the entries in the cluster)
# Tue  3 Feb 15:56:33 GMT 2026 1.7 added --topproteomes

# imports
import os
import sys
import time
import argparse
from tqdm import tqdm
from typing import List, Optional
from collections import Counter
from statistics import median


# constants
HEADER = "cluster_id\tprotein_ids\tproteins_count\tproteomes_count\trepresentative\tseqlen_range\tseqlen_mode\n"
DESCRIPTION = """
Script to extract lists of protein identifiers for clusters defined in a tsv file
    (like the output files of the label_clusters script)
The file is assumed to contain at least three columns, with the first three
being "cluster_id", "protein_id" and "proteomes".

This script will check the clusters in the input file, optionally filter only those
containing a certain number of distinct proteomes, and write the list of proteins
belonging to the cluster, each tagged with the proteome the protein comes from (only
the first one if many are listed in the proteomes column unless --all is passed)
into a space separated string.

The last three columns of the output will contain the first protein_id of each cluster (assumed to
be the representative one), the range of sequence lengths observed (if provided in protein_id in
the form protein_id|seqlen; the first entry will be checked to determine this) and the sequence length mode (the most frequent sequence length in the cluster).

Note: if there is no unique mode (if more than one length shares the maximum frequency), the rounded median is instead used (this is indicated by using a float instead of an integer).

Sample input file:
    cluster_id\tprotein_id\tproteomes\tis_rep
    0\tENSSSCP00000055324|661\t35497\t*
    0\tENSSSCP00055011301|568\t4698922\t
    0\tENSSSCP00035021320|596\t4698918\t

Sample output file:
    cluster_id\tprotein_ids\tproteins_count\tproteomes_count\trepresentative\tseqlen_range\tseqlen_mode
    0\t35497:ENSSSCP00000055324|661 4698918:ENSSSCP00035021320|596 4698922:ENSSSCP00055011301|568\t3\t3\tENSSSCP00000055324|661\t568-661\t661

Example call:
    ./filter_clusters.py --input_file results_pig/Labelled_Specie_protein_cluster.tsv --out_file results_pig/Protein_clusters_m13.tsv --minproteomes 13 --all -q
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

    parser = CustomArgumentParser(description="Protein clusters filter.")
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
        "-m",
        "--minproteomes",
        type=positive_integer,
        required=False,
        default=1,
        help="Optionally filter clusters which contain proteins from at least this number of unique proteomes",
    )
    parser.add_argument(
        "-t",
        "--topproteomes",
        type=positive_integer,
        required=False,
        help="Optionally filter clusters which contain proteins from no more than this number of unique proteomes",
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
        "-s",
        "--strict",
        action="store_true",
        required=False,
        dest="strict",
        help="Ignore any protein_id which is not labelled as belonging to a proteome_id",
        default=False,
    )
    parser.add_argument(
        "-a",
        "--all",
        action="store_true",
        required=False,
        dest="allproteomes",
        help="Keep all the proteome tags if several are listed in the proteomes column of the input file. Otherwise only the first one will be printed",
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
    eprint(f" |-- min proteomes threshold: {args.minproteomes}")
    if args.topproteomes is not None:
        eprint(f" |-- top proteomes threshold: {args.topproteomes}")

    return args


def filter_clusters(
    input_file,
    out_file,
    min_threshold=1,
    top_threshold=None,
    allproteomes=False,
    progress=False,
    strict=False,
):
    """
    Parse clusters in the input_file and prints list of proteins belonging to each cluster on
    a separate line, tagging them with the proteome(s) they belong to (only the first one or
    all of them according to allproteomes flag).

    Optionally, only those clusters which have at least a certain number of proteomes in it will be
    printed.

    Args:
        input_file (str): Path to the labelled clusters file.
        out_file (str): Path to the combined output file.

    Returns:
        int: The total number of proteins and clusters written to the output file.
    """
    protein_count = -1  # number of total proteins
    cluster_count = -1  # number of total clusters
    printed_protein_count = 0  # number of total proteins in the printed clusters
    printed_cluster_count = 0  # number of total clusters printed
    this_cluster_id = -1  # cluster_id of the cluster we are processing
    this_cluster_proteins = []  # proteins in the cluster we are processing
    this_line_proteomes = []  # proteomes found in each line
    this_cluster_proteomes = set()  # unique set of proteomes in the current cluster
    this_cluster_proteins_count = 0  # number of proteins in the current cluster
    this_cluster_lengths = []  # sequence lengths in the cluster
    proteomes_per_cluster_sum = 0  # to compute average of proteomes per cluster
    representative = ""  # stores the representative (first) protein_id of the cluster
    min_seq_len = 1e10
    max_seq_len = -1e10
    seqlen_range = ""
    seqlen_mode = ""
    possible_seqlen = False

    input_file_size = os.path.getsize(input_file)
    if progress:
        pbar = tqdm(total=input_file_size, unit="B")

    with open(out_file, "w") as output_fh:
        output_fh.write(HEADER)
        with open(input_file, "r") as input_fh:
            # read first two lines separately
            # in order to be able to skip header and get first representative efficiently
            line1 = input_fh.readline()
            line2 = input_fh.readline()
            columns = line2.rstrip("\n").split("\t", maxsplit=3)
            cluster_id, protein_id, proteome_ids = columns[0:3]
            representative = protein_id
            if "|" in protein_id:
                parts = protein_id.split("|")
                if len(parts) > 1:
                    try:
                        seqlen = int(parts[1])
                        # mark that we can probably get seqlen information after identifier
                        possible_seqlen = True
                        if seqlen < min_seq_len:
                            min_seq_len = seqlen
                        if seqlen > max_seq_len:
                            max_seq_len = seqlen
                        this_cluster_lengths.append(seqlen)
                    except:
                        pass

            # return to beginning of file
            input_fh.seek(0)
            if line1.startswith("cluster_id"):
                next(input_fh)  # skip the header only if present

            for line in input_fh:
                protein_count += 1
                if progress:
                    pbar.update(len(line))
                columns = line.rstrip("\n").split("\t", maxsplit=3)
                # if len(columns) < 3:
                #   eprint(f"found line with less than 3 columns: {line}") #debug
                #   continue
                cluster_id, protein_id, proteome_ids = columns[0:3]
                if strict and proteome_ids == "":
                    eprint(
                        f"skipping {protein_id} from {cluster_id} because it has no proteome label (strict mode)"
                    )
                    continue
                this_line_proteomes = proteome_ids.split(",")  # could be several
                if not allproteomes:
                    this_line_proteomes = this_line_proteomes[
                        0:1
                    ]  # keep only first one if several
                if (
                    possible_seqlen
                    and cluster_id == this_cluster_id
                    and "|" in protein_id
                ):
                    parts = protein_id.split("|")
                    if len(parts) > 1:
                        try:
                            seqlen = int(parts[1])
                            if seqlen < min_seq_len:
                                min_seq_len = seqlen
                            if seqlen > max_seq_len:
                                max_seq_len = seqlen
                            this_cluster_lengths.append(seqlen)
                        except:
                            pass
                if cluster_id != this_cluster_id:
                    cluster_count += 1
                    # we went through one full cluster, let's print it out
                    this_cluster_proteomes_count = len(this_cluster_proteomes)
                    this_cluster_proteins_count = len(this_cluster_proteins)
                    proteomes_per_cluster_sum += this_cluster_proteomes_count

                    if min_seq_len != 1e10:  # if we got some seqlen information
                        seqlen_range = f"{min_seq_len}-{max_seq_len}"
                        seqlen_counts = Counter(this_cluster_lengths)
                        seqlen_modes = [
                            v
                            for v, f in seqlen_counts.items()
                            if f == max(seqlen_counts.values())
                        ]
                        seqlen_mode = (
                            seqlen_modes[0]
                            if len(seqlen_modes) == 1
                            else f"{median(this_cluster_lengths):.1f}"
                        )

                    if this_cluster_proteomes_count >= min_threshold and (
                        top_threshold is None
                        or this_cluster_proteomes_count <= top_threshold
                    ):
                        output_fh.write(
                            f"{this_cluster_id}\t{' '.join(sorted(this_cluster_proteins))}\t{this_cluster_proteins_count}\t{len(this_cluster_proteomes)}\t{representative}\t{seqlen_range}\t{seqlen_mode}\n"
                        )
                        printed_protein_count += this_cluster_proteins_count
                        printed_cluster_count += 1
                    # else: eprint(f"skipping printing {this_cluster_id} because it contains proteins from only {this_cluster_proteomes_count} proteomes: {this_cluster_proteomes}") #debug

                    # start a new cluster with the current line
                    this_cluster_id = cluster_id
                    representative = protein_id  # representative for this new cluster
                    this_cluster_proteins = set(
                        [
                            proteome_id + ":" + protein_id
                            for proteome_id in this_line_proteomes
                        ]
                    )
                    this_cluster_proteomes = set(this_line_proteomes)
                    # reset:
                    min_seq_len = 1e10
                    max_seq_len = -1e10
                    this_cluster_proteins_count = 0
                    this_cluster_lengths = []
                    seqlen_range = ""
                    seqlen_mode = ""
                    if possible_seqlen and "|" in protein_id:
                        parts = protein_id.split("|")
                        if len(parts) > 1:
                            try:
                                seqlen = int(parts[1])
                                if seqlen < min_seq_len:
                                    min_seq_len = seqlen
                                if seqlen > max_seq_len:
                                    max_seq_len = seqlen
                                this_cluster_lengths.append(seqlen)
                            except:
                                pass
                    # eprint(f"starting a new cluster {cluster_id} with {this_cluster_proteins}") #debug
                    continue
                # else: still in the same cluster
                # eprint(f"adding {[proteome_id + ':' + protein_id for proteome_id in this_line_proteomes]} to cluster {cluster_id}") #debug
                this_cluster_proteins.update(
                    [
                        proteome_id + ":" + protein_id
                        for proteome_id in this_line_proteomes
                    ]
                )
                this_cluster_proteomes.update(this_line_proteomes)

            # print the last cluster
            protein_count += 1
            cluster_count += 1
            # eprint(f"print final cluster {cluster_id} with {this_cluster_proteins_count} proteins") #debug
            this_cluster_proteomes_count = len(this_cluster_proteomes)
            this_cluster_proteins_count = len(this_cluster_proteins)
            proteomes_per_cluster_sum += this_cluster_proteomes_count

            if min_seq_len != 1e10:  # if we got some seqlen information
                seqlen_range = f"{min_seq_len}-{max_seq_len}"
                seqlen_counts = Counter(this_cluster_lengths)
                seqlen_modes = [
                    v
                    for v, f in seqlen_counts.items()
                    if f == max(seqlen_counts.values())
                ]
                seqlen_mode = (
                    seqlen_modes[0]
                    if len(seqlen_modes) == 1
                    else f"{median(this_cluster_lengths):.1f}"
                )

            if this_cluster_proteomes_count >= min_threshold:
                output_fh.write(
                    f"{this_cluster_id}\t{' '.join(sorted(this_cluster_proteins))}\t{this_cluster_proteins_count}\t{len(this_cluster_proteomes)}\t{representative}\t{seqlen_range}\t{seqlen_mode}\n"
                )
                printed_protein_count += this_cluster_proteins_count
                printed_cluster_count += 1

    if progress:
        pbar.close()
    return (
        protein_count,
        cluster_count,
        printed_protein_count,
        printed_cluster_count,
        proteomes_per_cluster_sum,
    )


if __name__ == "__main__":
    initial_secs = time.time()  # for total time count
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())

    eprint(f" .-- BEGUN {timestamp} --.")
    args = check_args()
    eprint(f" |...")

    (
        total_proteins,
        total_clusters,
        printed_proteins,
        printed_clusters,
        total_proteomes_per_cluster,
    ) = filter_clusters(
        args.input_file,
        args.out_file,
        min_threshold=args.minproteomes,
        top_threshold=args.topproteomes,
        allproteomes=args.allproteomes,
        progress=args.progress,
        strict=args.strict,
    )

    eprint(
        f" |-- seen {total_proteins} protein_ids in {total_clusters} clusters (mean of {total_proteins/total_clusters:.2f} proteins/cluster and {total_proteomes_per_cluster/total_clusters:.2f} proteomes/cluster)"
    )
    if printed_clusters:
        eprint(
            f" |-- printed {printed_proteins} protein_ids in {printed_clusters} clusters (mean of {printed_proteins/printed_clusters:.2f} proteins/cluster)"
        )
        eprint(
            " |-- final file created -- Elapsed: {}, {} clusters/s --".format(
                *elapsed_time(initial_secs, total_clusters),
            )
        )
    else:
        eprint(f" |-- nothing printed")
        delete_files([args.out_file])  # delete the file which has now only header
    # total time
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
    eprint(f" '-- ENDED {timestamp} --'")
