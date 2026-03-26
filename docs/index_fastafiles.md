
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

| Options | Values  | Help |
| ------- | ------- | ---- |
| <pre>-h --help</pre> | Flag. | show this help message and exit |
| <pre>-i --input\_file</pre> | Optional.<br/>Type: file | Path to the fasta file to be indexed |
| <pre>-l --list</pre> | Optional.<br/>Type: file | Path to a file containing a list of filepaths to index |
| <pre>-d --fasta\_dir</pre> | Optional.<br/>Type: str | Directory containing fasta files \(with .fa extension unless -e specified\) |
| <pre>-e --extension</pre> | Type: str<br/>Default: `.fa` | Extension for files in fasta\_dir\; default \'.fa\' |
| <pre>-f --force</pre> | Flag. | Force re-creation of existing index files |
| <pre>-q --progress</pre> | Flag. | Show a progress bar |
| <pre>-p --pattern</pre> | Default: `\^\(.\+\\\|.\*\)\$` | Regexp Capture pattern for identifiers to index in fasta header. Default is\: \'\^\(.\+\\\|.\*\)\$\' |
| <pre>-t --threads</pre> | Type: int<br/>Default: `1` | Number of threads for parallel processing |
