#!/usr/bin/env python3
# changelog
# Wed  4 Feb 10:26:15 GMT 2026 0.1 started
# Wed  4 Feb 11:20:55 GMT 2026 0.9 coded
# Wed  4 Feb 12:10:54 GMT 2026 1.0 working
# Wed  4 Feb 12:23:05 GMT 2026 1.1 compute_jaccard multi (not only best match)
# Wed  4 Feb 16:58:57 GMT 2026 1.2 added jaccardmin
# Wed  4 Feb 17:29:15 GMT 2026 1.3 improved formatting
# Thu  5 Feb 08:55:16 GMT 2026 1.4 compute_jaccard symmetric
# Wed 18 Feb 17:17:25 GMT 2026 1.5 analyze_events and write_events instead of analyze_splits_merges to deal also with complex cases
# Wed 18 Feb 21:36:08 GMT 2026 1.6 added --skip_metrics
# Tue 24 Feb 11:32:27 GMT 2026 1.7 added -swm -swo
# Fri 27 Feb 13:27:42 GMT 2026 1.8 correct printing of which outfiles get created
# Fri  6 Mar 11:08:56 GMT 2026 1.9 added IoU for proteins overlap
# Fri  6 Mar 11:33:11 GMT 2026 2.0 added containmentmin for Szymkiewicz–Simpson coefficient

# Script to compare two clusterings (e.g. obtained from different clustering parameters).
# 1) reads clusters from two Protein_clusters tsv files
#   creating cluster->protein and protein->cluster dictionaries
# 2) computes proteins overlap (like "comm": protein identifiers only present in file1,
#   common for the two files, only present in file2)
# 3) computes several sklearn clustering similarity metrics
# 4) uses Jaccard index (IoU) to create best mappings between clusters,
#   all those beyond a certain IoU threshold
# 4.2) uses Szymkiewicz–Simpson coefficient as alternative to Jaccard index to deal with asymmetric
#   sets, with overlap/containment threshold
# 5) checks the mappings to identify possible split or merge events
#   (single cluster of file1 whose proteins appear in multiple clusters of file2 and vice versa)
# 6) creates four output files for further analyses, for example with inspect_merge_splits.py


import os
import sys
import csv
import time
import argparse
import numpy as np
from collections import defaultdict, Counter, deque
from sklearn.metrics import (
    adjusted_rand_score,
    adjusted_mutual_info_score,
    # customise with other metrics as needed
)


TSV_HAS_HEADER = True
csv.field_size_limit(sys.maxsize)


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

    def zero_to_one(value):
        try:
            value = float(value)
            if value < 0:
                raise argparse.ArgumentTypeError(
                    "{} is not a number in interval [0.0-1.0]".format(value)
                )
            elif value > 1:
                raise argparse.ArgumentTypeError(
                    "{} is not a number in interval [0.0-1.0]".format(value)
                )
        except ValueError:
            raise Exception("{} is not a number in interval [0.0-1.0]".format(value))
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

    parser = argparse.ArgumentParser(description="Compare two clustering TSV files")

    parser.add_argument(
        "-m",
        "--minproteomes",
        type=positive_integer,
        default=1,
        help="Minimum number of proteomes per cluster",
    )

    parser.add_argument(
        "-j",
        "--jaccardmin",
        type=zero_to_one,
        default=0.3,
        help="Minimum Jaccard index to consider for mapping clusters, default 0.3",
    )

    parser.add_argument(
        "-c",
        "--containmentmin",
        type=zero_to_one,
        default=0.9,
        help="Minimum containment threshold (Szymkiewicz–Simpson coefficient) for mapping clusters due to overlap (good for asymmetric sized clusters): default 0.9",
    )

    parser.add_argument(
        "-t",
        "--topproteomes",
        type=positive_integer,
        help="Maximum number of proteomes per cluster",
    )

    parser.add_argument(
        "-p",
        "--pidtype",
        type=positive_integer,
        default=1,
        help="What protein identifier to consider, 1: full source_id|seqlen|upi information (default); 2: only source_id|seqlen; 3: only upi|seqlen (if upi available, else source_id)",
    )

    parser.add_argument(
        "file1", type=is_valid_file, help="First protein clusters TSV file"
    )

    parser.add_argument(
        "file2", type=is_valid_file, help="Second protein clusters TSV file"
    )

    parser.add_argument("-o", "--outprefix", default="clucmp", help="Output prefix")

    parser.add_argument(
        "-sm",
        "--skip_metrics",
        action="store_true",
        required=False,
        help="Skip similarity metrics computation",
        default=False,
    )
    parser.add_argument(
        "-swo",
        "--skip_writing_overlap",
        action="store_true",
        required=False,
        help="Skip writing overlap",
        default=False,
    )
    parser.add_argument(
        "-swm",
        "--skip_writing_mapping",
        action="store_true",
        required=False,
        help="Skip writing mapping",
        default=False,
    )

    args = parser.parse_args()
    eprint(f" |-- file1: '{args.file1}'")
    eprint(f" |-- file2: '{args.file2}'")
    eprint(f" |-- pidtype: {args.pidtype}")
    if args.minproteomes > 1:
        eprint(f" |-- min proteomes threshold: {args.minproteomes}")
    if args.topproteomes is not None:
        eprint(f" |-- top proteomes threshold: {args.topproteomes}")

    outfiles = " |-- outfiles:"
    if not args.skip_metrics:
        outfiles += f" {args.outprefix}.metrics.tsv"
    if not args.skip_writing_overlap:
        outfiles += f" {args.outprefix}.overlap.tsv"
    if not args.skip_writing_mapping:
        outfiles += f" {args.outprefix}.mapping.tsv"

    outfiles += f" {args.outprefix}.mergesplits.tsv"
    eprint(outfiles)
    eprint(" |--")

    return args


# data loading functions
def clean_pids(proteins, pidtype=1):
    """
    Remove proteome_id identifiers and clean pids,
    depending on what protein identifier to consider, 1: full source_id|seqlen|upi information (default); 2: only source_id; 3: only upi (if available, else source_id)
    Assumes protein identifiers in the form "SOURCE_ID|SEQLEN|UPI"
    Returns a set, so no duplicated identifiers
    """
    cleaned_proteins = []
    for protein in proteins:
        proteome, protein = protein.split(":")
        if pidtype != 1:  # for 1 do not modify pid
            parts = protein.split("|")
            if pidtype == 2:  # keep only source_id
                protein = "|".join(parts[0:2])
            elif pidtype == 3:
                if len(parts) > 2:  # keep only upi, fallback to pid
                    protein = f"{parts[2]}|{parts[1]}"  # upi
                else:
                    protein = "|".join(parts[0:2])

        cleaned_proteins.append(protein)

    return set(cleaned_proteins)  # squash together identical pids


def load_clusters(
    filename,
    min_threshold=1,
    top_threshold=None,
    pidtype=1,
):
    """
    Load input TSV file, filtering clusters if required,
    retaining pid information according to pidtype

    Returns:
        clusters: dict[cluster_id] -> set(proteins)
        protein_to_cluster: dict[protein] -> cluster_id
    """
    clusters = {}
    protein_to_cluster = {}

    with open(filename) as fh:
        reader = csv.DictReader(
            fh,
            delimiter="\t",
            fieldnames=[
                "cluster_id",
                "protein_ids",
                "proteins_count",
                "proteomes_count",
                "representative",
                "seqlen_range",
                "seqlen_mode",
            ],
        )
        if TSV_HAS_HEADER:
            next(reader, None)

        for row in reader:
            cid = row["cluster_id"]
            proteomes_count = int(row["proteomes_count"])

            if proteomes_count >= min_threshold and (
                top_threshold is None or proteomes_count <= top_threshold
            ):
                proteins = row["protein_ids"].split()
                proteins = clean_pids(proteins, pidtype=pidtype)
                clusters[cid] = proteins
                for p in proteins:
                    protein_to_cluster[p] = cid

    return clusters, protein_to_cluster


# protein overlap functions
def compute_protein_overlap(map1, map2):
    """
    Compute protein overlap information
    """

    p1 = set(map1.keys())
    p2 = set(map2.keys())

    only1 = p1 - p2
    only2 = p2 - p1
    common = p1 & p2

    return only1, common, only2


def write_protein_overlap(only1, common, only2, outprefix, map1, map2):
    """
    Write protein overlap information
    """

    out = f"{outprefix}.overlap.tsv"

    with open(out, "w") as fh:
        fh.write("category\tprotein_id\tcluster_ids\n")
        for p in sorted(only1):
            fh.write(f"only_file1\t{p}\t{map1[p]}\n")
        for p in sorted(common):
            fh.write(f"common\t{p}\t{map1[p]},{map2[p]}\n")
        for p in sorted(only2):
            fh.write(f"only_file2\t{p}\t{map2[p]}\n")

    # eprint(f" |->   Protein overlap written to {out}")


# similarity metrics functions
def build_label_vectors(map1, map2, common):
    """
    Build aligned label arrays for sklearn metrics
    """
    labels1 = []
    labels2 = []

    for p in common:
        labels1.append(int(map1[p]))
        labels2.append(int(map2[p]))

    return labels1, labels2


def compute_similarity(labels1, labels2):
    """
    Compute a series of clustering similarity metrics
    """
    ari = adjusted_rand_score(labels1, labels2)
    ami = adjusted_mutual_info_score(labels1, labels2)

    return {
        "ARI": ari,
        "AMI": ami,
    }


def write_similarity(scores, outprefix):
    """
    Write the similarity scores
    """
    out = f"{outprefix}.metrics.tsv"
    with open(out, "w") as fh:
        fh.write("metric\tvalue\n")
        for k, v in scores.items():
            fh.write(f"{k}\t{v:.6f}\n")

    # eprint(f" |->   Similarity scores written to {out}")


# Jaccard mapping functions
def compute_jaccard_mapping_symmetric(
    clusters1,
    clusters2,
    map1,
    map2,
    common,
    min_jaccard=0.2,
    min_containment=0.8,  # containment threshold for Szymkiewicz–Simpson coefficient
):
    """
    Compute a symmetric Jaccard-based mapping between two clusterings.

    This function builds the full (sparse) contingency matrix between
    clusters in file1 and file2 by iterating over shared proteins and
    counting co-memberships.

    For each protein present in both clusterings:
        - Determine its cluster in file1 and in file2
        - Increment the intersection count for that cluster pair

    This yields all non-zero intersections (|A ∩ B|) between clusters.

    For each cluster pair (A, B) with non-zero overlap:
        - Compute the Jaccard index:
              J(A,B) = |A ∩ B| / (|A| + |B| - |A ∩ B|)
        - Retain pairs whose Jaccard value meets the minimum threshold

    This approach is symmetric with respect to file order and produces
    a bipartite overlap graph between clusterings, suitable for detecting
    both split and merge events.

    Clusters that share no proteins are not represented, as their Jaccard
    index is zero by definition.

    Returns all significant cluster overlaps.
    """
    intersections = Counter()
    # Build contingency matrix
    for p in common:
        c1 = map1[p]
        c2 = map2[p]
        intersections[(c1, c2)] += 1

    mapping = []

    for (c1, c2), inter in intersections.items():
        size1 = len(clusters1[c1])
        size2 = len(clusters2[c2])
        union = size1 + size2 - inter
        jaccard = inter / union
        containment = inter / (min(size1, size2))

        if jaccard >= min_jaccard or containment >= min_containment:
            mapping.append(
                (
                    c1,
                    c2,
                    size1,
                    size2,
                    inter,
                    union,
                    jaccard,
                    containment,
                )
            )

    return mapping


def write_mapping(mapping, outprefix):
    """
    Write jaccard mapping
    """
    out = f"{outprefix}.mapping.tsv"

    with open(out, "w") as fh:
        fh.write("clu1\tclu2\tsize1\tsize2\tI\tU\tIoU\toverl\n")
        for row in mapping:
            fh.write("{}\t{}\t{}\t{}\t{}\t{}\t{:.6f}\t{:.4f}\n".format(*row))
    # eprint(f" |->   Jaccard mapping written to {out}")


def analyze_events(mapping):
    """
    Analyze merge/split/complex events using connected components
    of the bipartite cluster overlap graph.
    e.g.
    A B → X (merge)
    A → X Y (split)
    A B → X Y (complex)

    Returns:
        events: list of (c1_list, c2_list)
        counts: [n_splits, n_merges, n_complex]
    """
    graph = defaultdict(set)

    # build undirected bipartite graph
    for c1, c2, *_ in mapping:
        graph[("f1", c1)].add(("f2", c2))
        graph[("f2", c2)].add(("f1", c1))

    seen = set()
    merges = []
    splits = []
    complex_events = []

    for node in graph:
        if node in seen:
            continue

        queue = deque([node])
        seen.add(node)

        c1s = set()
        c2s = set()

        while queue:
            side, cid = queue.popleft()

            if side == "f1":
                c1s.add(cid)
            else:
                c2s.add(cid)

            for nb in graph[(side, cid)]:
                if nb not in seen:
                    seen.add(nb)
                    queue.append(nb)

        n1 = len(c1s)
        n2 = len(c2s)

        # skip one-to-one
        if n1 == 1 and n2 == 1:
            continue

        c1s = sorted(c1s)
        c2s = sorted(c2s)

        if n1 > 1 and n2 == 1:
            merges.append((c1s, c2s))
        elif n1 == 1 and n2 > 1:
            splits.append((c1s, c2s))
        else:
            complex_events.append((c1s, c2s))

    # final ordered list
    events = merges + splits + complex_events

    counts = [
        len(splits),
        len(merges),
        len(complex_events),
    ]

    return events, counts


def write_events(events, outprefix):
    """
    Write merge/split/complex events to file
    """

    out = f"{outprefix}.mergesplits.tsv"

    with open(out, "w") as fh:
        fh.write("file1\tfile2\n")
        for c1s, c2s in events:
            fh.write(f"{' '.join(c1s)}\t{' '.join(c2s)}\n")


# Main
def compare_clusters(args):
    # note clusters are uniqued, without repeated identifiers
    clusters1, map1 = load_clusters(
        args.file1,
        min_threshold=args.minproteomes,
        top_threshold=args.topproteomes,
        pidtype=args.pidtype,
    )

    clusters2, map2 = load_clusters(
        args.file2,
        min_threshold=args.minproteomes,
        top_threshold=args.topproteomes,
        pidtype=args.pidtype,
    )

    eprint(f" |-- Clusters file1: {len(clusters1)}")
    eprint(f" |-- Clusters file2: {len(clusters2)}")

    # Protein overlap
    only1, common, only2 = compute_protein_overlap(map1, map2)
    only1_count = len(only1)
    common_count = len(common)
    only2_count = len(only2)
    eprint(f" |-- Proteins overlap")
    eprint(f" |--   only_file1: {only1_count}")
    eprint(f" |--   common:     {common_count}")
    eprint(f" |--   only_file2: {only2_count}")
    eprint(f" |--   IoU: {common_count/(only1_count + common_count + only2_count):.2f}")
    if not args.skip_writing_overlap:
        write_protein_overlap(only1, common, only2, args.outprefix, map1, map2)

    if not common:
        exit_with_error("[ERROR] No common proteins. Exiting.", 1)

    # Rest of the analysis is only on common set

    # Jaccard/Szymkiewicz–Simpson mapping
    mapping = compute_jaccard_mapping_symmetric(
        clusters1,
        clusters2,
        map1,
        map2,
        common,
        min_jaccard=args.jaccardmin,
        min_containment=args.containmentmin,
    )
    eprint(f" |-- Mapped clusters: {len(mapping)}")
    eprint(f" |-- (Jaccard min threshold: {args.jaccardmin})")
    eprint(f" |-- (Containment min threshold: {args.containmentmin})")
    if not args.skip_writing_mapping:
        write_mapping(mapping, args.outprefix)

    # Clusters merge/split events
    # splits, merges = analyze_splits_merges(mapping)
    # eprint(f" |-- Splits: {len(splits)}")
    # eprint(f" |-- Merges: {len(merges)}")
    # write_splits_merges(splits, merges, args.outprefix)
    events, counts = analyze_events(mapping)
    eprint(f" |-- Splits: {counts[0]}")
    eprint(f" |-- Merges: {counts[1]}")
    eprint(f" |-- Complex: {counts[2]}")
    write_events(events, args.outprefix)

    if not args.skip_metrics:
        # Similarity metrics
        labels1, labels2 = build_label_vectors(map1, map2, common)
        scores = compute_similarity(labels1, labels2)
        eprint(f" |-- Similarity metrics:")
        eprint(f" |--   {({k: round(v, 3) for k, v in scores.items()})}")
        write_similarity(scores, args.outprefix)


if __name__ == "__main__":
    initial_secs = time.time()  # for total time count
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())

    eprint(f" .-- BEGUN {timestamp} --.")

    args = check_args()

    compare_clusters(args)
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
    eprint(f" '-- ENDED {timestamp} --'")
