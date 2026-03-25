#!/usr/bin/env python3
# changelog
# Tue 17 Feb 12:08:34 GMT 2026 0.1 started
# Tue 17 Feb 15:37:55 GMT 2026 1.0 coded first version
# Wed 18 Feb 12:24:58 GMT 2026 1.1 made more efficient using only 4 subprocess (parallel) calls
# Wed 18 Feb 16:37:19 GMT 2026 1.2 added printing of mapping information with cluster sizes and IoU after parsing mapping file for relevant lines

# This script can be used to inspect cluster merge/split events, from the output of compare_clusters.py
# it parses the merge/split file
# optionally subsamples groups (--maxgroups)
# for each selected group:
# - extracts clusters from file1 and file2
# - runs alignments on the extracted FASTA files

# For performance reasons, we will use only 4 subprocess calls:
# Across all selected groups:
# Collect all unique cluster IDs needed from column1 (file1 clusters) and column2 (file2 clusterds)
# Then run extract step once per file and align step once per file
# Finally loop over groups to print recap paths

import os
import sys
import random
import argparse
import subprocess
from multiprocessing import cpu_count
from collections import defaultdict
from pathlib import Path


UPIAPIURL = (
    "https://rest.uniprot.org/uniparc/{UPI}.fasta"  # to retrieve sequences by UPI
)
BASEPATH = "/homes/insana/proteomescomparisons/github/"  # TMP
EXTRACT_CLUSTERS = f"{BASEPATH}scripts/extract_clusters.py"
ALIGN_CLUSTERS = f"{BASEPATH}scripts/dedup_align_fastafiles.py"


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


def run_cmd(cmd):
    # print("[CMD]", " ".join(cmd))
    # subprocess.run(cmd, check=True, stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)
    subprocess.run(cmd, check=True)


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

    parser = argparse.ArgumentParser(
        description="Inspect cluster merge/split events by extracting clusters and aligning them"
    )

    parser.add_argument(
        "mergesplits",
        type=is_valid_file,
        help="Merge/split file from compare_clusters.py (two tab-separated columns, space-separated cluster IDs), e.g. '90vs80.mergesplits.tsv'",
    )
    parser.add_argument(
        "mapping",
        type=is_valid_file,
        help="Cluster mapping TSV from compare_clusters.py, e.g. '90vs80.mapping.tsv'",
    )
    parser.add_argument(
        "file1",
        type=is_valid_file,
        help="First clustering TSV file, e.g. '90/Protein_clusters_m.tsv'",
    )
    parser.add_argument(
        "file2",
        type=is_valid_file,
        help="Second clustering TSV file, e.g. '80/Protein_clusters_m.tsv'",
    )
    parser.add_argument(
        "-os",
        "--outsuffix",
        default="clusters",
        help="Suffix for output directories (default: clusters)",
    )
    parser.add_argument(
        "-op",
        "--outprefix",
        default="",
        help="Optional prefix for output directories, e.g. 'dataset'",
    )
    parser.add_argument(
        "-m",
        "--maxsample",
        type=positive_integer,
        help="Randomly analyze up to this many merge/split groups",
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=positive_integer,
        required=False,
        default=1,
        help="Number of threads for parallel processing",
    )

    return parser.parse_args()


def load_mapping(mapping_path, file1_clusters):
    """
    Load cluster mapping file and index by clu1
    Load only mapping rows with clu1 in the provided list
    Returns: dict[clu1] -> list of rows
    """

    mapping = defaultdict(list)
    with open(mapping_path) as fh:
        header = fh.readline().strip().split("\t")
        for line in fh:
            if not line.strip():
                continue

            fields = line.rstrip().split("\t")

            clu1 = fields[0]
            if clu1 in file1_clusters:
                mapping[clu1].append(line)

    return mapping


def inspect_groups(args):
    """
    for each selected group:
        collect unique cluster IDs needed from column1 (file1 clusters) and column2 (file2 clusterds)
    then run extract step once per file and align step once per file
    finally loop over groups to print recap paths
    """
    mergesplits_path = Path(args.mergesplits)
    file1 = Path(args.file1)
    file2 = Path(args.file2)

    groups = []
    file1_cluster_ids = set()
    file2_cluster_ids = set()

    with open(mergesplits_path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith("file1") or line.startswith("#") or not line:
                continue

            parts = line.split("\t")
            if len(parts) != 2:
                exit_with_error(
                    f"ERROR: malformed line found in mergesplits input: {line}", 5
                )
                continue

            c1 = parts[0].split()
            c2 = parts[1].split()

            groups.append((c1, c2))

    if not groups:
        exit_with_error("ERROR: No valid merge/split groups found", 34)

    if args.threads > cpu_count():
        args.threads = cpu_count()
        eprint(f" |-- WARNING: only {args.threads} threads available")
    if args.threads > 1:
        eprint(f" |-- threads: {args.threads}")
    args.threads = str(args.threads)

    # subsample if requested
    if args.maxsample and len(groups) > args.maxsample:
        eprint(
            f" |-- Subsampling {args.maxsample} out of {len(groups)} merge/split groups"
        )
        groups = random.sample(groups, args.maxsample)
    else:
        eprint(f" |-- Processing all {len(groups)} merge/split groups")

    for c1, c2 in groups:
        file1_cluster_ids.update(c1)
        file2_cluster_ids.update(c2)

    # output directories
    prefix = f"{args.outprefix}_" if args.outprefix else ""

    file1_out = Path(f"{prefix}file1{args.outsuffix}")
    file2_out = Path(f"{prefix}file2{args.outsuffix}")

    file1_out.mkdir(parents=True, exist_ok=True)
    file2_out.mkdir(parents=True, exist_ok=True)

    # extract all clusters
    eprint(f" |-- Extracting {len(file1_cluster_ids)} clusters from file1")
    run_cmd(
        [
            EXTRACT_CLUSTERS,
            "-i",
            str(file1),
            "-o",
            str(file1_out),
            "-c",
            " ".join(file1_cluster_ids),
            "-w",
            UPIAPIURL,
            "-q",
            "-ut",  # only unique sequences and terse header
        ]
    )
    eprint(f" |-- Extracting {len(file2_cluster_ids)} clusters from file2")
    run_cmd(
        [
            EXTRACT_CLUSTERS,
            "-i",
            str(file2),
            "-o",
            str(file2_out),
            "-c",
            " ".join(file2_cluster_ids),
            "-w",
            UPIAPIURL,
            "-q",
            "-ut",  # only unique sequences and terse header
        ]
    )

    # align all clusters
    eprint(" |-- Aligning all file1 clusters")
    run_cmd(
        [
            ALIGN_CLUSTERS,
            "-o",
            str(file1_out),
            "--outfmt",
            "clustal",
            "-d",
            str(file1_out),
            "-t",
            args.threads,
            "-a",  # skip deduplication
            "-q",
        ]
    )
    eprint(" |-- Aligning all file2 clusters")
    run_cmd(
        [
            ALIGN_CLUSTERS,
            "-o",
            str(file2_out),
            "--outfmt",
            "clustal",
            "-d",
            str(file2_out),
            "-t",
            args.threads,
            "-a",  # skip deduplication
            "-q",
        ]
    )

    # load mapping information
    mapping_by_clu1 = load_mapping(Path(args.mapping), file1_cluster_ids)

    # final recap
    eprint("\n |-- MERGE / SPLIT INSPECTION RECAP --|")

    for c1, c2 in groups:
        is_merge = len(c1) > 1 and len(c2) == 1
        is_split = len(c1) == 1 and len(c2) > 1

        if is_merge:
            eventlabel = "Merge"
        elif is_split:
            eventlabel = "Split"
        else:
            eventlabel = "Complex"

        # check if we have all aln files for them
        aln_files = []

        for cid in c1:
            aln = file1_out / f"{cid}.aln"
            if aln.exists():
                aln_files.append(str(aln))

        for cid in c2:
            aln = file2_out / f"{cid}.aln"
            if aln.exists():
                aln_files.append(str(aln))

        if len(aln_files) == len(c1) + len(c2):
            lhs = " ".join(c1)
            rhs = " ".join(c2)
            print(f"\n{eventlabel} {lhs} -> {rhs}")
            print("clu1\tclu2\tsize1\tsize2\tI\tU\tIoU")

            # print mapping info
            for clu1 in c1:
                for row in mapping_by_clu1.get(clu1, []):
                    print(row, end="")

            print("    " + " ".join(aln_files))
        # otherwise we skip, as at least one failed alignment (probably single unique sequence)


if __name__ == "__main__":
    args = check_args()
    inspect_groups(args)
