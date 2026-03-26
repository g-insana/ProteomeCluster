
This script extracts FASTA sequences for protein identifiers grouped by cluster from a TSV (like the output of filter_clusters.py), retrieving sequences from ffdb-indexed FASTA files and writing one FASTA file per cluster.

 * Expects a TSV with at least two columns: cluster_id and protein_ids
   - Each entry in the protein_ids column is a space-separated token formatted as: proteome_id:protein_id
 * Writes, for each cluster, a FASTA file named <cluster_id>.fa containing the sequences of the proteins in that cluster
   - Sequences are extracted from FASTA files located in FASTA_DIR using ffdb index files (filename.idx)
   - Proteome IDs map to FASTA filenames; the script supports recovering filenames using optional prefixes/extensions


Sample input file:
    cluster_id	protein_ids	proteins_count	proteomes_count
    0	35497:ENSSSCP00000055324|661 4698922:ENSSSCP00055011301|568 4698918:ENSSSCP00035021320|596	3	3

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

usage: extract_clusters.py [-h] (-d FASTA_DIR | -s SINGLE_FILE | -w URL) -i INPUT_FILE -o OUT_DIR [-u]
                           [-a] [-ut] [-m MINPROTEINS] [-r RESTRICT] [-c CLUSTER [CLUSTER ...]]
                           [-p PREFIX] [-e EXTENSION] [-t THREADS] [-f] [-q]

Clusters sequence extraction.

optional arguments:
  -h, --help            show this help message and exit
  -p PREFIX, --prefix PREFIX
                        Prefix for fasta filenames, to be added to proteome identifiers; e.g.
                        'proteome_'
  -e EXTENSION, --extension EXTENSION
                        Extension for files in fasta_dir, to be added to proteome identifiers; e.g.
                        '.fa'
  -t THREADS, --threads THREADS
                        Number of threads for parallel processing.
  -f, --force           Force overwrite of existing cluster output files.
  -q, --progress        Show a progress bar

Fasta sequences from (choose one):
  -d FASTA_DIR, --fasta_dir FASTA_DIR
                        Directory containing all proteome fasta files
  -s SINGLE_FILE, --single_file SINGLE_FILE
                        Path to a single fasta file where to extract all identifiers from
  -w URL, --web URL     URL to web API where to fetch fasta sequences by UPI identifier e.g.
                        'https://rest.uniprot.org/uniparc/{UPI}.fasta' ({UPI} is where the UniParc
                        identifier will be placed in the request

Input and Output:
  -i INPUT_FILE, --input_file INPUT_FILE
                        Path to the input file
  -o OUT_DIR, --out_dir OUT_DIR
                        Directory where the cluster .fa files will be written

Comprehensive output:
  -u, --uniq            Only keep unique sequences (merging the identifiers into one header) in case
                        multiple protein identifiers have the same sequence.
  -a, --all             Keep all proteome labels in case the same protein_id is tagged to multiple
                        proteomes in the input file.
  -ut, --uniqterse      Only keep unique sequences (under the first seen identifier and without
                        proteome label: shortest header)

Filtering options:
  -m MINPROTEINS, --minproteins MINPROTEINS
                        Optionally only extract clusters that contain at least this number of
                        proteins. Note that this is checked first, before applying restricted filter.
  -r RESTRICT, --restrict RESTRICT
                        Path to a file containing proteome identifiers. If given, only sequences
                        belonging to the proteomes from that file will be extracted.
  -c CLUSTER [CLUSTER ...], --cluster CLUSTER [CLUSTER ...]
                        Optionally specify which cluster(s) to extract, by cluster_id (multiple
                        cluster_id can be specified, space separated)
