
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
    cluster_id	protein_id	proteomes	is_rep
    0	ENSSSCP00000055324|661	35497	*
    0	ENSSSCP00055011301|568	4698922	
    0	ENSSSCP00035021320|596	4698918	

Sample output file:
    cluster_id	protein_ids	proteins_count	proteomes_count	representative	seqlen_range	seqlen_mode
    0	35497:ENSSSCP00000055324|661 4698918:ENSSSCP00035021320|596 4698922:ENSSSCP00055011301|568	3	3	ENSSSCP00000055324|661	568-661	661

Example call:
    ./filter_clusters.py --input_file results_pig/Labelled_Specie_protein_cluster.tsv --out_file results_pig/Protein_clusters_m13.tsv --minproteomes 13 --all -q

usage: filter_clusters.py [-h] -i INPUT_FILE -o OUT_FILE [-m MINPROTEOMES] [-t TOPPROTEOMES] [-q] [-s]
                          [-a]

Protein clusters filter.

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_FILE, --input_file INPUT_FILE
                        Path to the input file
  -o OUT_FILE, --out_file OUT_FILE
                        Path to the output file
  -m MINPROTEOMES, --minproteomes MINPROTEOMES
                        Optionally filter clusters which contain proteins from at least this number of
                        unique proteomes
  -t TOPPROTEOMES, --topproteomes TOPPROTEOMES
                        Optionally filter clusters which contain proteins from no more than this
                        number of unique proteomes
  -q, --progress        Show a progress bar
  -s, --strict          Ignore any protein_id which is not labelled as belonging to a proteome_id
  -a, --all             Keep all the proteome tags if several are listed in the proteomes column of
                        the input file. Otherwise only the first one will be printed
