
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

usage: index_fastafiles.py [-h] (-i INPUT_FILE | -l LIST | -d FASTA_DIR) [-e EXTENSION] [-f] [-q]
                           [-p PATTERN] [-t THREADS]

Indexer for fasta files.

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_FILE, --input_file INPUT_FILE
                        Path to the fasta file to be indexed
  -l LIST, --list LIST  Path to a file containing a list of filepaths to index
  -d FASTA_DIR, --fasta_dir FASTA_DIR
                        Directory containing fasta files (with .fa extension unless -e specified)
  -e EXTENSION, --extension EXTENSION
                        Extension for files in fasta_dir; default '.fa'
  -f, --force           Force re-creation of existing index files
  -q, --progress        Show a progress bar
  -p PATTERN, --pattern PATTERN
                        Regexp Capture pattern for identifiers to index in fasta header. Default is:
                        '^(.+\|.*)$'
  -t THREADS, --threads THREADS
                        Number of threads for parallel processing
