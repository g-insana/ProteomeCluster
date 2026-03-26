
This script labels the second column of a TSV (protein IDs) with the FASTA filename(s) (proteome identifiers) where each protein appears, compacts cluster IDs to sequential integers, and flags the representative protein for each cluster. Input row order is preserved.

 * Expects a two-column TSV (like mmseqs2 cluster files): column 1 = cluster identifier, column 2 = protein identifier (matching FASTA headers).
 * Scans FASTA files in FASTA_DIR and maps protein IDs (from headers) to the FASTA filename
 * Adds a new column, proteomes, containing the proteome identifier(s) for each protein. If a protein is found in multiple FASTA files, identifiers are comma-separated.
 * Replaces original cluster identifiers with sequential integers starting at 0 while preserving the original row order and grouping identical original cluster IDs to the same integer.
 * Adds a fourth column, is_rep, marking the first protein encountered for each original cluster with '*' (representative); other rows in that cluster have an empty value.
 * Produces an output TSV with header:
    cluster_id	protein_id	proteomes	is_rep

Required arguments:
    FASTA_DIR  - directory containing FASTA files to search (headers provide protein IDs)
    INPUT_FILE - input TSV file (two columns: cluster identifier and protein identifier)
    OUT_FILE   - output TSV to create (four columns)

Sample input file (cluster and protein columns (tab-separated):
    ENSSSCP00000055324|661	ENSSSCP00000055324|661
    ENSSSCP00000055324|661	ENSSSCP00055011301|568
    ENSSSCP00000055324|661	ENSSSCP00035021320|596

Sample output file (header included):
    cluster_id	protein_id	proteomes	is_rep
    0	ENSSSCP00000055324|661	35497	*
    0	ENSSSCP00055011301|568	4698922	
    0	ENSSSCP00035021320|596	4698918	

Example call:
    label_clusters.py --fasta_dir pig/ --input_file results_pig/Specie_protein_cluster.tsv --out_file results_pig/Labelled_Specie_protein_cluster.tsv --prefix proteome_ --extension .fa

Distributed approach (batch tagging, optionally parallel):
      When dealing with a huge amount of entries memory requirements could be a problem.
      The labelling work could then be performed in batches, specifying a --batchsize. That will be the
      maximum number of fasta_files to read in one go to create a partial dictionary and annotate the
      input file with that. The input file will be annotated several times until all fasta_files have
      been ingested.

      This can be combined with --threads to do the tagging in parallel.
      The input file will be split into thrice the number of threads and each thread will only read
      a number of fasta files equal to batchsize at each time.

      Optionally --chunksize can be specified: the input file will be split in chunks of that size
          (minimum 5Mb)

      E.g.: ./label_clusters.py --fasta_dir ecoli/ --input_file results_ecoli/Specie_protein_cluster.tsv --out_file results_ecoli/Labelled_Specie_protein_cluster.tsv --prefix proteome_ --extension .fa --threads 10 --batchsize 1000 --chunksize 1G

      Alternative approach: the input file will be copied and processed independently by each worker and
      results combined in the end. To use this approach specify --chunksize=n.
      This approach eliminates the chance of time spent with workers trying to acquire file lock on
      the same chunk, but it uses more disk space and could suffer i/o performance loss when combining
      result files in the end.

usage: label_clusters.py [-h] -f FASTA_DIR -i INPUT_FILE -o OUT_FILE [-p PREFIX] [-e EXTENSION]
                         [-n NOLABEL] [-s] [-u] [-t THREADS] [-b BATCHSIZE] [-c CHUNKSIZE] [-q]

Proteome labeller for clusters of protein identifiers.

optional arguments:
  -h, --help            show this help message and exit
  -f FASTA_DIR, --fasta_dir FASTA_DIR
                        Directory containing all .fa files
  -i INPUT_FILE, --input_file INPUT_FILE
                        Path to the input file
  -o OUT_FILE, --out_file OUT_FILE
                        Path to the output file
  -p PREFIX, --prefix PREFIX
                        Optionally prefix for filenames, which will be removed; e.g. 'proteome_'
  -e EXTENSION, --extension EXTENSION
                        Extension for files in fasta_dir (e.g. '.fa')
  -n NOLABEL, --nolabel NOLABEL
                        Optional string for missing mapping; default is '' (e.g. '?')
  -s, --sortlabels      Optionally sort the labels in the output file
  -u, --uniq            Optionally uniq the input file (in case it has repeated lines)
  -t THREADS, --threads THREADS
                        Number of threads for parallel processing
  -b BATCHSIZE, --batchsize BATCHSIZE
                        Max number of fasta files that will be processed at the same time
  -c CHUNKSIZE, --chunksize CHUNKSIZE
                        Chunk size in which to split the input file; if not specified: split the input
                        file in a number of chunks equal to thrice the number of threads; minimum
                        chunksize: 5m; use 'n' to avoid splitting the input file
  -q, --progress        Show a progress bar
