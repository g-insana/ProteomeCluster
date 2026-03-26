
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

usage: dedup_align_fastafiles.py [-h]
                                 (-d INPUT_FASTA_DIR | -l INPUT_FASTA_LIST | -s FASTA_FILE [FASTA_FILE ...])
                                 [-t THREADS] [-e EXTENSION] [-a] [-q] [-m] [-c] -o OUT_DIR
                                 [-at ALIGN_THREADS] [-f] [-st {Protein,DNA}]
                                 [-of {fasta,clustal,stockholm}]

Clustal Omega wrapper with sequence deduplication

optional arguments:
  -h, --help            show this help message and exit
  -t THREADS, --threads THREADS
                        Number of threads for parallel processing.
  -e EXTENSION, --extension EXTENSION
                        File extension for FASTA files (default: .fa).
  -a, --all             Skip the deduplication step and align all sequences
  -q, --progress        Show a progress bar
  -m, --mapfiles        Create mapping files with duplicated sequences
  -c, --consensus       Writes consensus line to file .cns, if present (clustal format only)

Input Options (choose one):
  -d INPUT_FASTA_DIR, --fasta_dir INPUT_FASTA_DIR
                        Directory containing fasta files to be aligned.
  -l INPUT_FASTA_LIST, --fasta_list INPUT_FASTA_LIST
                        Path to a text file containing a list of FASTA files.
  -s FASTA_FILE [FASTA_FILE ...], --single_file FASTA_FILE [FASTA_FILE ...]
                        Fasta file name(s), space separated

Clustal Omega specific arguments:
  -o OUT_DIR, --out_dir OUT_DIR
                        Output directory for aligned files.
  -at ALIGN_THREADS, --align_threads ALIGN_THREADS
                        Number of threads for Clustal Omega.
  -f, --force           Force overwrite of existing output files.
  -st {Protein,DNA}, --seqtype {Protein,DNA}
                        Specify the sequence type (e.g., 'Protein' or 'DNA').
  -of {fasta,clustal,stockholm}, --outfmt {fasta,clustal,stockholm}
                        Specify the MSA output format.
