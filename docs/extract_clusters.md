
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

| Options | Values  | Help |
| ------- | ------- | ---- |
| <pre>-h --help</pre> | Flag. | show this help message and exit |
| <pre>-p --prefix</pre> | Optional.<br/>Type: str | Prefix for fasta filenames, to be added to proteome identifiers\; e.g. \'proteome\_\' |
| <pre>-e --extension</pre> | Optional.<br/>Type: str | Extension for files in fasta\_dir, to be added to proteome identifiers\; e.g. \'.fa\' |
| <pre>-t --threads</pre> | Default: `1` | Number of threads for parallel processing. |
| <pre>-f --force</pre> | Flag. | Force overwrite of existing cluster output files. |
| <pre>-q --progress</pre> | Flag. | Show a progress bar |
| *Fasta sequences from \(choose one\)* | |
| <pre>-d --fasta\_dir</pre> | Optional.<br/>Type: str | Directory containing all proteome fasta files |
| <pre>-s --single\_file</pre> | Optional.<br/>Type: file | Path to a single fasta file where to extract all identifiers from |
| <pre>-w --web</pre> | Optional.<br/>Type: str | URL to web API where to fetch fasta sequences by UPI identifier e.g. \'https\://rest.uniprot.org/uniparc/\{UPI\}.fasta\' \(\{UPI\} is where the UniParc identifier will be placed in the request |
| *Input and Output* | |
| <pre>-i --input\_file</pre> | INPUT_FILE<br/>Required.<br/>Type: file | Path to the input file |
| <pre>-o --out\_dir</pre> | OUT_DIR<br/>Required.<br/>Type: str | Directory where the cluster .fa files will be written |
| *Comprehensive output* | |
| <pre>-u --uniq</pre> | Flag. | Only keep unique sequences \(merging the identifiers into one header\) in case multiple protein identifiers have the same sequence. |
| <pre>-a --all</pre> | Flag. | Keep all proteome labels in case the same protein\_id is tagged to multiple proteomes in the input file. |
| <pre>-ut --uniqterse</pre> | Flag. | Only keep unique sequences \(under the first seen identifier and without proteome label\: shortest header\) |
| *Filtering options* | |
| <pre>-m --minproteins</pre> | Type: int<br/>Default: `1` | Optionally only extract clusters that contain at least this number of proteins. Note that this is checked first, before applying --restrict filter. |
| <pre>-r --restrict</pre> | Optional.<br/>Type: file | Path to a file containing proteome identifiers. If given, only sequences belonging to the proteomes from that file will be extracted. |
| <pre>-c --cluster</pre> | Optional.<br/>Type: str | Optionally specify which cluster\(s\) to extract, by cluster\_id \(multiple cluster\_id can be specified, space separated\) |
