# ProtClustComp - Proteome Clustering & Comparisons Framework

Pipelines and tools for clustering, comparing and analyzing proteomes at scale.

ProtClustComp combines a **Nextflow pipeline** with **Python utilities** to perform clustering, filtering, alignments and other comparative analyses of proteomes.

---

## 🚀 Overview

ProtClustComp enables:

- High-performance **clustering** using MMseqs2
- Flexible **cluster filtering and annotation**
- Optional **sequence extraction and multiple alignment**
- Advanced **comparison of clustering strategies**
- Exploration of **merge/split cluster events**
- Analyses of **protein frequency** across isolates
- Inter- and intra- cluster **sequence variation**
- **Variome analysis**

The pipeline is designed to run both **locally** and on **HPC systems (SLURM/LSF)**.
The single processes can also be run independently via the provided python scripts.

---

## ✨ Features

### 🔹 Basic Workflow

- Clusters proteins using `mmseqs2`
- Labels cluster members by proteome of origin
- Filters clusters based on proteome representation

**Outputs:**
(for each species or dataset)
- `*_clusters.tsv`
  - Cluster membership
  - Sequence length range and mode
- `*_matrix.tsv`
  - presence absence matrix with row and column counts

---

### 🔹 Extended Workflow

Optional downstream processing after clustering:

- Extract sequences per cluster
- Align them by cluster

**Outputs:**

- Cluster FASTA files
- Cluster alignment files

---

### 🔹 Meta-Clustering Comparison

Compare different clustering runs (e.g. parameter tuning):

- Compute similarity metrics
- Identify overlapping clusters
- Detect cluster **merge/split events**

---

### 🔹 Variation Analysis (WIP)

- Framework for building and analyzing **variomes**

---

## 📦 Requirements

- [MMseqs2](https://github.com/soedinglab/MMseqs2)
- [Nextflow](https://www.nextflow.io/)
- Python ≥ 3.8

Optional:

- [ffdb](https://github.com/g-insana/ffdb.py) (for sequence extraction)
- [Biopython](https://biopython.org/) (for alignment processing)
- Clustal Omega (`clustalo`) (for multiple sequence alignments)
- [scikit-learn](https://scikit-learn.org) (for meta-clustering metrics)

---

## 📥 Installation

```bash
git clone https://github.com/g-insana/ProtClustComp.git
```

---

## ⚡ Quickstart 1 2 3

### 1. Prepare Input

Place proteome FASTA files under:
```
proteomesdir/
```

You can organize them in subfolders (e.g. per species or per dataset) and they will be processed in parallel.

---

### 2. Run Pipeline
either locally or on HPC, choosing your preferred parameters, editing the following scripts - or [running directly from command line](#command-line)

#### Basic workflow

```bash
./cluster_only-local.sh
./cluster_only-slurm.sh
```

#### Extended workflow

```bash
./cluster_and_align-local.sh
./cluster_and_align-slurm.sh
```

---

### 3. Results

Output files will be written to:
```
outdir/
```

---

### 🧪 Test Dataset

Included:

- 10 proteomes from [UniProtKB](https://https://www.uniprot.org)
- belonging to two species

Run immediately without setup.

---

## Command Line

Instead of the provided scripts you can for example run it on command line, e.g.:

```bash
nextflow run cluster_and_align.nf \
  --min_seq_id 0.9 \
  --coverage 0.9 \
  --proteomes_threshold 3 \
  -profile local -resume
```

See full list of parameters below in the [Documentation](#-documentation) section.

---

## 🧠 Best Practices

- Start with **strict thresholds** (for coverage, sequence identity, covmode) for high-confidence clusters
- Use **meta-clustering tools** (see below for the Documentation) to tune parameters and compare with the high-confidence clustering
- Align `maxForks` and threads on HPC systems to your resources and needs
- Use `--extract` and `--align` only if needed

---

## 📚 Documentation

### Nextflow Pipeline
#### Core Parameters

| Parameter | Description |
|----------|------------|
| `--proteomesdir` | Input proteome directory |
| `--outdir` | Output directory |
| `--proteomes_threshold` | Minimum proteomes per cluster |
| `--extract` | Extract sequences |
| `--align` | Perform alignments |

---

#### MMseqs2 Clustering

| Parameter | Description |
|----------|------------|
| `--covmode` | Coverage mode (0=bidirectional, 1=target, 2=query) |
| `--coverage` | Minimum alignment coverage (0–1) |
| `--min_seq_id` | Minimum sequence identity (0–1) |
| `--seq_id_mode` | Identity calculation mode |

---

#### Performance

| Parameter | Description |
|----------|------------|
| `--clusterthreads` | Threads for clustering |
| `--labelthreads` | Threads for labeling |
| `--extractthreads` | Threads for extraction |
| `--alignthreads` | Threads for alignment |
| `--clustalothreads` | Threads for Clustal Omega |
| `--maxForks` | Max parallel processes |

---

#### Advanced

```bash
--mmseqs_additional "--clust-hash 1 --kmer-per-seq 100 ..."
```
(custom MMseqs2 parameters for tuning speed vs sensitivity vs accuracy)

### Cluster Comparison Tools

Two scripts are provided:

### 1️⃣ compare_clusters.py

Compare two clustering results:
- compute overlap between protein identifiers of the two clusterings
- map clusters according to Jaccard and Szymkiewicz–Simpson thresholds (IoU and containment)
- write the mapping between clusters of the two clustering results
- identify split/merge/complex events
- compute similarity metrics (`adjusted_rand_score` and `adjusted_mutual_info_score`)

```bash
compare_clusters.py clusters_file1.tsv clusters_file2.tsv
```

#### Key Options

| Option | Description |
|-------|------------|
| `-m` | Minimum proteomes per cluster |
| `-j` | Minimum Jaccard index (default 0.3) |
| `-c` | Containment threshold: Szymkiewicz–Simpson coefficient (default 0.9) |
| `-t` | Max proteomes per cluster |
| `-p` | Protein ID format |
| `-o` | Output prefix |

#### Optional Flags

- `--skip_metrics`
- `--skip_writing_overlap`
- `--skip_writing_mapping`

---

### 2️⃣ inspect_merge_splits.py

Analyze cluster merge/split events:
- extract clusters involved (sampling if desired)
- align them
- prepare list of alignments for visual comparison

```bash
inspect_merge_splits.py mergesplits.tsv mapping.tsv file1.tsv file2.tsv
```

#### Key Options (all optional)

| Option | Description |
|-------|------------|
| `-os` | Output suffix |
| `-op` | Output prefix |
| `-m` | Max sampled events |
| `-t` | Threads |


### Variation Analysis
(WIP)

### Individual processes

Instead of running the Nextflow pipeline you can use the individual Python scripts to execute specific steps independently. This is useful for debugging, custom workflows, or integrating ProtClustComp components into other pipelines.

Below is a brief documentation for each script. Complete documentation is available (simply click on the script name below).

---

## 🏷️ Labelling

### [`label_clusters.py`](docs/label_clusters.md)

This script labels the second column of a TSV (protein IDs) with the FASTA filename(s) (proteome identifiers) where each protein appears, compacts cluster IDs to sequential integers, and flags the representative protein for each cluster. Input row order is preserved.

```bash
usage: label_clusters.py [-h] -f FASTA_DIR -i INPUT_FILE -o OUT_FILE [-p PREFIX] [-e EXTENSION]
                         [-n NOLABEL] [-s] [-u] [-t THREADS] [-b BATCHSIZE] [-c CHUNKSIZE] [-q]
```

---

## 🔍 Filtering

### [`filter_clusters.py`](docs/filter_clusters.md)

Script to extract lists of protein identifiers for clusters defined in a tsv file
    (like the output files of `label_clusters.py`)
The file is assumed to contain at least three columns, with the first three
being `cluster_id`, `protein_id` and `proteomes`.


```bash
usage: filter_clusters.py [-h] -i INPUT_FILE -o OUT_FILE [-m MINPROTEOMES] [-t TOPPROTEOMES] [-q] [-s]
                          [-a]
```

---

## 📦 Extraction

### [`extract_clusters.py`](docs/extract_clusters.md)

This script extracts FASTA sequences for protein identifiers grouped by cluster from a TSV (like the output of `filter_clusters.py`), retrieving sequences from ffdb-indexed FASTA files and writing one FASTA file per cluster.

```bash
usage: extract_clusters.py [-h] (-d FASTA_DIR | -s SINGLE_FILE | -w URL) -i INPUT_FILE -o OUT_DIR [-u]
                           [-a] [-ut] [-m MINPROTEINS] [-r RESTRICT] [-c CLUSTER [CLUSTER ...]]
                           [-p PREFIX] [-e EXTENSION] [-t THREADS] [-f] [-q]

```

---

## 🧬 Alignment

### [`dedup_align_fastafiles.py`](docs/dedup_align_fastafiles.md)

Wrapper around Clustal Omega that can find and deduplicate sequences across FASTA files, run alignments and optionally create mapping files of duplicates.

```bash
usage: dedup_align_fastafiles.py [-h]
                                 (-d INPUT_FASTA_DIR | -l INPUT_FASTA_LIST | -s FASTA_FILE [FASTA_FILE ...])
                                 [-t THREADS] [-e EXTENSION] [-a] [-q] [-m] [-c] -o OUT_DIR
                                 [-at ALIGN_THREADS] [-f] [-st {Protein,DNA}]
                                 [-of {fasta,clustal,stockholm}]
```

---

## 🧮 Presence/Absence Matrix

### [`presence_matrix.py`](docs/presence_matrix.md)

Script to create a presence/absence (and count) matrix where the columns
are proteomes, the rows are clusters and the values are either 0 if there are no proteins from a
column proteome in the row cluster or a positive integer corresponding to the number of proteins from
that proteome in the cluster.

```bash
usage: presence_matrix.py [-h] -i INPUT_FILE -p PROTEOMES -o OUT_FILE [-m MEMBERS_FILE] [-l LABELS]
                          [-a ADDITIONAL] [-t] [-c] [-s]
```

---

## ⚡FASTA Indexing (Extra Utility)

### [`index_fastafiles.py`](docs/index_fastafiles.md)

Indexes FASTA files for fast random access. It supports indexing a single file, a list of files or all the files with a given extension present in a directory.

```bash
usage: index_fastafiles.py [-h] (-i INPUT_FILE | -l LIST | -d FASTA_DIR) [-e EXTENSION] [-f] [-q]
                           [-p PATTERN] [-t THREADS]
```

---

## 📜 License

ProtClustComp is licensed under the [Apache License 2.0](https://www.apache.org/licenses/LICENSE-2.0.html)

(c) Copyright [Giuseppe Insana](https://insana.net), 2024-

---

## 🤝 Contributions

Contributions, issues, and feature requests are welcome.

---

## 📌 Citation

If you use ProtClustComp in your research, please consider citing this repository.

---
