
Wrapper around Clustal Omega that can find and deduplicate sequences across FASTA files, run alignments and optionally create mapping files of duplicates.

 * Accepts either a directory of FASTA files or a text file listing FASTA paths (mutually exclusive; one required)
 * Deduplicates identical sequences across input FASTA files (optionally recording mapping files of duplicated sequence occurrences (-m))
 * Alternatively align all sequences (-a), skipping deduplication step
 * Runs Clustal Omega to produce multiple sequence alignments for each input FASTA, using a configurable number of threads
 * Writes aligned output files into a specified output directory (created if missing)
 * Tweakable parallelism: a global worker pool (--threads) manages file-level parallelism; Clustal Omega itself receives a separate thread count (--align_threads) per job
 * Optionally write the Clustal consensus line in a separate file for each alignment

Example calls:
    dedup_align_fastafiles.py -d fasta_dir -o aligned_out -t 4 -at 2 # align all FASTAs in a directory with 4 workers and 2 threads per Clustal Omega run
    dedup_align_fastafiles.py -l my_fastas.txt -o aligned_out -m -q # align files listed in a text file, create duplicate mapping files, show progress
    dedup_align_fastafiles.py -l my_fastas.txt -a -o aligned_out # align all sequences of listed files, skipping deduplication step

| Options | Values  | Help |
| ------- | ------- | ---- |
| *positional arguments* | |
| *optional arguments* | |
| <pre>-h --help</pre> | Flag. | show this help message and exit |
| <pre>-t --threads</pre> | Type: Unknown<br/>Default: `1` | Number of threads for parallel processing. |
| <pre>-e --extension</pre> | Type: str<br/>Default: `.fa` | File extension for FASTA files \(default\: .fa\). |
| <pre>-a --all</pre> | Flag. | Skip the deduplication step and align all sequences |
| <pre>-q --progress</pre> | Flag. | Show a progress bar |
| <pre>-m --mapfiles</pre> | Flag. | Create mapping files with duplicated sequences |
| <pre>-c --consensus</pre> | Flag. | Writes consensus line to file .cns, if present \(clustal format only\) |
| *Input Options \(choose one\)* | |
| <pre>-d --fasta\_dir</pre> | Optional.<br/>Type: str | Directory containing fasta files to be aligned. |
| <pre>-l --fasta\_list</pre> | Optional.<br/>Type: Unknown | Path to a text file containing a list of FASTA files. |
| <pre>-s --single\_file</pre> | Optional.<br/>Type: str | Fasta file name\(s\), space separated |
| *Clustal Omega specific arguments* | |
| <pre>-o --out\_dir</pre> | OUT_DIR<br/>Required.<br/>Type: str | Output directory for aligned files. |
| <pre>-at --align\_threads</pre> | Type: Unknown<br/>Default: `1` | Number of threads for Clustal Omega. |
| <pre>-f --force</pre> | Flag. | Force overwrite of existing output files. |
| <pre>-st --seqtype</pre> | Type: str<br/>Choice: `Protein`, `DNA`<br/>Default: `Protein` | Specify the sequence type \(e.g., \'Protein\' or \'DNA\'\). |
| <pre>-of --outfmt</pre> | Type: str<br/>Choice: `fasta`, `clustal`, `stockholm`<br/>Default: `fasta` | Specify the MSA output format. |
