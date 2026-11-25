#!/usr/bin/env python3
# Sun  7 Dec 14:58:12 GMT 2025 0.1
# Sun  7 Dec 21:59:53 GMT 2025 0.9
# Mon  8 Dec 21:42:11 GMT 2025 1.0
# Tue  9 Dec 12:17:52 GMT 2025 1.1
# Thu  8 Jan 09:58:23 GMT 2026 1.2 added --strict
# Thu  8 Jan 11:26:20 GMT 2026 1.3 added --counts
# Wed 14 Jan 13:31:13 GMT 2026 1.4 added merging of clusters having same label
# Fri 13 Feb 12:39:12 GMT 2026 1.5 added --members_file

import os
import sys
import time
import argparse
from collections import defaultdict


DESCRIPTION = """
Script to create a presence/absence (and count) matrix where the columns
are proteomes, the rows are clusters and the values are either 0 if there are no proteins from a
column proteome in the row cluster or a positive integer corresponding to the number of proteins from
that proteome in the cluster. The format will be saved in tsv format.

Minimal required input is a filtered protein clusters file (-i)
  (like the output of the filter_clusters script) whose first two columns contains cluster_id and protein_ids (space separated entries in the form proteome_id:protein_id)
and a file containing a list of proteome_ids (-p)
  (assumed to be all the proteome_ids present in the protein clusters file).

Sample input file:
    cluster_id\tprotein_ids\tproteins_count\tproteomes_count\trepresentative\tseqlen_range
    0\t35497:ENSSSCP00000055324|661 4698918:ENSSSCP00035021320|596 4698922:ENSSSCP00055011301|568\t3\t3\tENSSSCP00000055324|661\t568-661

Sample output file:
    cluster\t35497\t4229143\t4698918\t4698922
          0\t1\t0\t1\t1

If the proteomes list file contains two columns, the second column will be assumed to contain proteome labels which will be used to replace the proteome_ids from the output file as column labels.
The order of proteomes in the proteomes list file also determines the orders of the columns in the output matrix.

An optional clusterslabel file can be specified (-l), which will be used to name the rows. This file should contain two columns: the first one containing cluster_ids and the second one containing the cluster labels.

An optional additional file (-a) can be specified. This file is assumed to contain singleton clusters (clusters with protein(s) from a single proteome). The second column of this file will be considered to contain cluster_ids, the third column protein_id and the fourth column containing proteome_ids.

If the -t option is passed, a column "proteomes_count" will be added (as 2nd column) with the number of proteomes having at least one protein for the cluster corresponding to each row (same as the 4th column of the clusters input file).

If the -c option is passed, a final row "clusters_per_proteome" will be added with the number of clusters present for each proteome (number of clusters in which the proteome is present with at least one protein).

An optional cluster members file can be created, specifying its filename (-m).

Example call:
    presence_matrix.py -i Protein_clusters_m2.tsv -o matrix.tsv -p proteomelabels.tsv -l clusterlabels.tsv -a add_singletons.acc -t
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
        description="Create presence/absence TSV matrix from a filtered protein clusters tsv. Optionally create cluster members file as well."
    )
    parser.add_argument(
        "-i",
        "--input_file",
        type=is_valid_file,
        required=True,
        help="Path to the input file"
        "first column cluster_id, "
        "second column protein_ids (proteome_id:protein_id space separated pairs)",
    )
    parser.add_argument(
        "-p",
        "--proteomes",
        type=is_valid_file,
        required=True,
        help="Path to the list of proteomes or proteome labels mapping file: either one (proteome_ids) or two columns (proteome_id proteome_label)",
    )
    parser.add_argument(
        "-o",
        "--out_file",
        type=str,
        required=True,
        help="Path to the matrix output file",
    )
    parser.add_argument(
        "-m",
        "--members_file",
        type=str,
        required=False,
        help="Path to the optional cluster members file",
    )
    parser.add_argument(
        "-l",
        "--labels",
        required=False,
        help="optional cluster labels mapping file: two columns (cluster_id cluster_label); NOTE that if same label is applied to different cluster_id, the corresponding rows will be merged",
        type=is_valid_file,
    )
    parser.add_argument(
        "-a",
        "--additional",
        required=False,
        type=is_valid_file,
        help="optional .acc file for additional singleton clusters to add to the output matrix: 2nd col = cluster_id, 4th col = proteome_id",
    )
    parser.add_argument(
        "-t",
        "--totals",
        action="store_true",
        required=False,
        help="Add a column with total proteomes per cluster",
        default=False,
    )
    parser.add_argument(
        "-c",
        "--counts",
        action="store_true",
        required=False,
        help="Add a final row with cluster counts per proteome",
        default=False,
    )
    parser.add_argument(
        "-s",
        "--strict",
        action="store_true",
        required=False,
        help="Ensure totals match (if -t used)",
        default=False,
    )

    args = parser.parse_args()

    eprint(f" |-- input_file: {args.input_file}")

    try:
        with open(args.out_file, "w"):
            pass
        os.remove(args.out_file)
    except permissionerror:
        exit_with_error(f"error: cannot write to file '{args.out_file}'.", 1)
    eprint(f" |-- out_file: {args.out_file}")

    if args.members_file is not None:
        try:
            with open(args.members_file, "w"):
                pass
            os.remove(args.members_file)
        except permissionerror:
            exit_with_error(f"error: cannot write to file '{args.members_file}'.", 1)
        eprint(f" |-- members_file: {args.members_file}")

    if args.strict and not args.totals:
        exit_with_error(f"ERROR: --strict only allowed with --totals", 2)

    eprint(f" |-- proteomes list: {args.proteomes}")
    if args.labels is not None:
        eprint(f" |-- cluster labels from: {args.labels}")
    if args.additional is not None:
        eprint(f" |-- additional singleton clusters from: {args.additional}")
    return args


def load_proteome_labels(path):
    """
    proteome_id -> proteome_id (if a single column file is provided)

    or

    proteome_id -> proteome_label
    e.g. 18\tUP000000625
    """
    d = {}
    proteomes_order = []
    with open(path) as f:
        line1 = f.readline()  # header?
        f.seek(0)  # return to beginning of file
        if any(c.isalpha() for c in line1.split("\t", 1)[0]):
            eprint(" |-- NOTICE: skipping header in proteomes' labels file")
            next(f)  # skip the header if present
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")[:2]
            if len(parts) == 2:
                prot, label = parts
            else:
                prot = label = parts[0]
            d[prot] = label
            proteomes_order.append(label)
    return d, proteomes_order


def load_cluster_labels(path):
    """
    cluster_id -> cluster_label
    e.g. 77437\tUPI000012E897
    """
    mapping = {}  # cid_to_label
    label_to_cids = {}
    with open(path) as f:
        line1 = f.readline()  # header?
        f.seek(0)  # return to beginning of file
        if any(c.isalpha() for c in line1.split("\t", 1)[0]):
            eprint(" |-- NOTICE: skipping header in clusters' labels file")
            next(f)  # skip the header if present
        for line in f:
            if not line.strip():
                continue
            cid, label = line.rstrip("\n").split("\t")[:2]
            mapping[cid] = label
            if label not in label_to_cids:
                label_to_cids[label] = []
            label_to_cids[label].append(cid)

    reused_labels = {l: cids for l, cids in label_to_cids.items() if len(cids) > 1}

    return mapping, reused_labels


def parse_clusters(clusters_file, proteome_labels, cluster_labels, membersfh=None):
    """
    Parse the clusters file and yield:
        (cluster_label, dict of proteomes -> count)
    """
    with open(clusters_file) as f:
        line1 = f.readline()  # header?
        f.seek(0)  # return to beginning of file
        if line1.startswith("cluster_id"):
            eprint(" |-- NOTICE: skipping header in clusters file")
            next(f)  # skip the header if present
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                eprint(" => WARNING, wrong format for line in clusters file!")
            cid, proteins = parts[0:2]
            proteomes_count = parts[3]

            if cluster_labels == {}:
                cluster_label = cid  # use cluster_id if cluster labels not provided
            else:
                if cid not in cluster_labels:
                    # skip unknown cluster_ids
                    continue
                cluster_label = cluster_labels[cid]

            counts = defaultdict(int)

            if proteins:
                for entry in proteins.split():
                    if ":" not in entry:
                        continue
                    prot_id, protein_id = entry.split(":", 1)
                    if prot_id in proteome_labels:
                        proteome_label = proteome_labels[prot_id]
                        counts[proteome_label] += 1
                        if membersfh is not None:
                            membersfh.write(
                                f"{proteomes_count}\t{cluster_label}\t{proteome_label}\t{protein_id}\n"
                            )

            yield cid, cluster_label, counts, proteomes_count


def parse_additional(
    additional_file, proteome_labels, cluster_labels, seen_clusters, membersfh=None
):
    """
    Parse .acc file for singleton clusters not present in the .tsv file.
    Yields (cid, source_id, counts)
    """
    store = defaultdict(lambda: defaultdict(int))

    with open(additional_file) as f:
        line1 = f.readline()  # header?
        f.seek(0)  # return to beginning of file
        if any(c.isalpha() for c in line1.split("\t", 1)[0]):
            eprint(" |-- NOTICE: skipping header in additional_file")
            next(f)  # skip the header if present
        for line in f:
            if not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 4:
                eprint(" => WARNING, wrong format for line in additional file")
                continue

            cid, protein_id, proteome_id = cols[1:4]

            if cid in seen_clusters:
                eprint(
                    f" => NOTICE: ignoring cluster {cid} from additional file, already parsed from main file"
                )
                # already output from main .tsv file
                continue

            if cluster_labels == {}:
                cluster_label = cid
            else:
                if cid not in cluster_labels:
                    eprint(
                        f" => NOTICE: ignoring cluster {cid} from additional file, missing cluster_label"
                    )
                    # skip unknown cluster_ids
                    continue
                else:
                    cluster_label = cluster_labels[cid]

            if proteome_id in proteome_labels:
                proteome_label = proteome_labels[proteome_id]
                store[cid][proteome_label] += 1
                if membersfh is not None:
                    membersfh.write(
                        f"1\t{cluster_label}\t{proteome_label}\t{protein_id}\n"
                    )  # proteomes_count assumed to be 1 for additional
            else:
                eprint(
                    f" => NOTICE: ignoring cluster {cid} from additional file, unknown proteome identifier {proteome_id}"
                )

    for cid, proteome_counts in store.items():
        if cluster_labels == {}:
            cluster_label = cid
        else:
            cluster_label = cluster_labels[cid]
        yield cid, cluster_label, proteome_counts


if __name__ == "__main__":
    initial_secs = time.time()  # for total time count

    # timestamp = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
    # eprint(f" .-- BEGUN {timestamp} --.")
    args = check_args()

    # Load labels
    proteome_labels, all_proteome_labels = load_proteome_labels(args.proteomes)
    if args.labels is not None:
        cluster_labels, reused_labels = load_cluster_labels(args.labels)
    else:
        cluster_labels = reused_labels = {}

    any_reused_labels = len(reused_labels) > 0
    if any_reused_labels:
        eprint(
            f" |-- NOTICE: Reused cluster labels (will be joined): {list(reused_labels.keys())}"
        )

    merged_counts = {label: defaultdict(int) for label in reused_labels}

    membersfh = None
    if args.members_file is not None:
        membersfh = open(args.members_file, "w")
        membersfh.write("proteomes_count\tcluster\tproteome\tmember\n")

    with open(args.out_file, "w") as outfh:
        # Header
        if args.totals:
            outfh.write(
                "cluster\tproteomes_count\t" + "\t".join(all_proteome_labels) + "\n"
            )
        else:
            outfh.write("cluster\t" + "\t".join(all_proteome_labels) + "\n")

        if args.counts:
            clusters_per_proteome = [0] * len(all_proteome_labels)
        # Stream rows from main tsv
        seen = set()
        for cid, cluster_label, counts, proteomes_count in parse_clusters(
            args.input_file, proteome_labels, cluster_labels, membersfh
        ):
            seen.add(cid)

            if any_reused_labels and cluster_label in reused_labels:
                # eprint(f"holding info for {cluster_label}") #debug
                # hold info, do not write out until end
                for prot, v in counts.items():
                    merged_counts[cluster_label][prot] += v
                continue

            row_vals = [counts.get(u, 0) for u in all_proteome_labels]
            if args.counts:
                for i, v in enumerate(row_vals):
                    if v > 0:
                        clusters_per_proteome[i] += 1
            if args.totals:
                num_present = sum(1 for v in row_vals if v > 0)
                if args.strict and num_present != int(proteomes_count):
                    eprint(f"clusterid {cid} ({cluster_label}): {row_vals}")
                    exit_with_error(
                        f"ERROR: miscount in proteomes_count! Expected: {proteomes_count}, found: {num_present}",
                        34,
                    )
                outfh.write(
                    "\t".join(
                        [cluster_label, proteomes_count] + list(map(str, row_vals))
                    )
                    + "\n"
                )
            else:
                outfh.write(
                    "\t".join([cluster_label] + list(map(str, row_vals))) + "\n"
                )

        # Handle .acc singletons (if provided)
        if args.additional:
            for cid, cluster_label, counts in parse_additional(
                args.additional, proteome_labels, cluster_labels, seen, membersfh
            ):
                if any_reused_labels and cluster_label in reused_labels:
                    # eprint(f"holding additional info for {cluster_label}") #debug
                    # hold info, do not write out until end
                    for prot, v in counts.items():
                        merged_counts[cluster_label][prot] += v
                    continue

                row_vals = [counts.get(u, 0) for u in all_proteome_labels]
                if args.counts:
                    for i, v in enumerate(row_vals):
                        if v > 0:
                            clusters_per_proteome[i] += 1
                if args.totals:
                    num_present = sum(1 for v in row_vals if v > 0)
                    if args.strict and num_present != 1:
                        exit_with_error(
                            f"ERROR: miscount in proteomes_count! Expected: 1, found: {num_present}",
                            34,
                        )
                    outfh.write(
                        "\t".join([cluster_label, "1"] + list(map(str, row_vals)))
                        + "\n"
                    )
                else:
                    outfh.write(
                        "\t".join([cluster_label] + list(map(str, row_vals))) + "\n"
                    )

        if any_reused_labels:
            # now deal with cluster_id with same labels -> merging their information
            for label, counts in merged_counts.items():
                # eprint(f"merging info for {label}") #debug
                row_vals = [counts.get(u, 0) for u in all_proteome_labels]

                if args.counts:
                    for i, v in enumerate(row_vals):
                        if v > 0:
                            clusters_per_proteome[i] += 1

                if args.totals:
                    proteomes_count = sum(1 for v in row_vals if v > 0)
                    outfh.write(
                        "\t".join(
                            [label, str(proteomes_count)] + list(map(str, row_vals))
                        )
                        + "\n"
                    )
                else:
                    outfh.write("\t".join([label] + list(map(str, row_vals))) + "\n")

        if args.counts:
            if args.totals:
                outfh.write(
                    "\t".join(
                        ["#clustercounts:", "-"] + list(map(str, clusters_per_proteome))
                    )
                    + "\n"
                )
            else:
                outfh.write(
                    "\t".join(
                        ["#clusters_per_proteome:"]
                        + list(map(str, clusters_per_proteome))
                    )
                    + "\n"
                )

    if args.members_file is not None:
        membersfh.close()

    # timestamp = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
    # eprint(f" '-- ENDED {timestamp} --'")
