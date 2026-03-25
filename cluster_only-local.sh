#!/usr/bin/env bash
#covmode 0: bidirectional
#min_seq_id 0.9: 90% minimum sequence identity
#coverage 0.9: alignment covers at least 90% of both target and query
#proteomes_threshold: minimum number of proteomes required to be present in clusters
#clusterthreads, labelthreads, extractthreads, alignthreads, clustalothreads: threads for various processes
#maxForks 1: maximum number of concurrent processes

#input and output paths:
proteomesdir="$(pwd)/proteomes/"
outdir="$(pwd)/outdir/"

nextflow run cluster_and_align.nf \
  --proteomesdir ${proteomesdir} \
  --outdir ${outdir} \
  --covmode 0 \
  --min_seq_id 0.9 \
  --coverage 0.9 \
  --proteomes_threshold 3 \
  --clusterthreads 1 --labelthreads 1 --extractthreads 1 --alignthreads 1 --clustalothreads 1 \
  --maxForks 1 \
  -profile local -resume

