# ProtComp - Proteome Comparison Framework: Pipeline and tools for clustering and comparing proteomes

A NextFlow pipeline and Python scripts for comprehensive analyses of proteomes data, providing:
- Clustering via MMseqs2
- Labelling, filtering and aggregation of clusters
- Alignments and Variation analyses
- Meta clustering comparisons

## Features

* Basic workflow:
The pipeline uses `mmseqs2` to cluster sequences of proteomes, with all parameters easily configurable.
It labels the cluster members according to the proteomes where each protein originates.
A desired level of filtering can also be applied, to keep only clusters with proteins from a specified number of proteomes.

It creates:
- Clusters TSV with full membership and additional information (like range of sequence length variation or sequence length mode of the cluster)
- Presence absence matrix

* Extended workflow:
The pipeline can also be directed to extract the sequences of the clusters and optionally align them, in order to study intra-cluster variation across isolates.

It creates:
- Cluster fasta files
- Cluster alignments

* Meta clustering
Also included are scripts to compare two different clusterings (for example obtained with different parameters) and compute similarity, identify overlap, find matching clusters and identify and explore split/merge events

* Variation analysis
Creation of variome (WIP)

The pipeline can run locally or on HPC via e.g. slurm/lsf.
The single processes can also be run independently via the provided python scripts.

## Requirements

Clustering requires [MMseqs2](https://github.com/soedinglab/MMseqs2), which must be installed and accessible via the command line.
For extraction of sequences: [ffdb](https://github.com/g-insana/ffdb.py).
For alignments, it depends on [biopython](https://biopython.org/) and [clustalo](https://en.wikipedia.org/wiki/Clustal).
And, only for computing meta-clustering similarity metrics: [scikit-learn](https://scikit-learn.org).

## Download and installation

``` bash
pip3 install git+https://github.com/g-insana/ProtComp.git
```
  (from [GitHub](https://github.com/g-insana/ProtComp/)).


## Quickstart

1. Place fasta files of proteomes to be analysed in `proteomesdir/` folder (several subfolders can be created for parallel processing, e.g. one per species)

2. Simply run one of the provided shell scripts for either the basic workflow:

``` bash
./cluster_only-local.sh
./cluster_only-slurm.sh #if on hpc
```

or the extended one:
``` bash
./cluster_and_align-local.sh
./cluster_and_align-slurm.sh #if on hpc
```

to start the NextFlow pipeline.

3. Results will appear under `outdir/` folder.

Note that some test files (a total of ten [UniProtKB](https://www.uniprot.org/) proteomes belonging to two species) are already provided, so you can try out the execution without any setup.

All parameters are configurable, simply edit or adapt the provided shell script or run directly from command line with desired options, e.g.:
```bash
nextflow run cluster_and_align.nf --proteomes_threshold 10 --coverage 0.8 --extract
```

## Documentation
(WIP)

### For the NextFlow pipeline

### For Comparison of Clusterings

### For Variation analysis

## Copyright

`ProtComp` is licensed under the [Apache License 2.0](https://www.apache.org/licenses/LICENSE-2.0.html).

(c) Copyright [Giuseppe Insana](https://insana.net), 2024-
