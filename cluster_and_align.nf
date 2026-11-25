/*
NextFlow example pipeline to cluster proteomes (located in subfolders) via mmseqs2, create representative fasta, extract cluster sequences and align them.

Steps are optional and controlled by command line parameters.

NOTE: you may need to adapt script to your fasta headers or use proteome files with fasta headers like
>IDENTIFIER|SEQLEN
or
>IDENTIFIER|SEQLEN|UPI
/*

/* Changelog:
Sun  9 Nov 21:28:47 GMT 2025 0.1 started, gather_fasta
Mon 10 Nov 08:52:21 GMT 2025 0.2 index_fasta
Sun 23 Nov 19:03:12 GMT 2025 0.3 cluster_proteomes, label_clusters
Wed 26 Nov 17:15:42 GMT 2025 0.4 filter_clusters, extract_clusters, align_clusters
Tue 24 Mar 21:56:08 GMT 2026 0.5 mmseqs_additional and dynamic steps choice
Tue 24 Mar 23:23:51 GMT 2026 0.6 create_matrix
Wed 25 Nov 16:42:09 GMT 2026 1.1 converted to strict nextflow syntax (nextflow>=26.04)
Wed 25 Nov 17:26:11 GMT 2026 1.2 added create_repr
*/

/* Resources parameters:

job requirements to work with a species of ~20k bacterial proteomes:
30G mem & 16 cpus for Labeller jobs
230Gb mem & 32 cpus for mmseqs jobs

*/

//for testing purposes (increase as needed):
params.clusterthreads = 2
params.labelthreads = 2
params.extractthreads = 2
params.alignthreads = 2
params.clustalothreads = 2

//number of concurrent parallel jobs
params.maxForks = 10

//input files
params.proteomesdir = "${projectDir}/proteomes/"
//output dir
params.outdir = "${projectDir}/outdir/"

//behaviour (steps to perform)
params.repr = false // if true, extract clusters' representative sequences
params.extract = false // if true, extract all clusters' sequences
params.align = false // if true, also align clusters' sequences

//mmseqs parameters
params.mmseqs_exe = "mmseqs" //specify appropriate path to mmseqs executable, if not in $PATH
params.covmode = 0 //0) bidirectional, 1) target coverage, 2) query coverage
params.seq_id_mode = 0 //0: alignment length 1: shorter, 2: longer sequence
params.coverage = 0.9 //[0.0-1.0] alignment covers at least __% of target (for cov-mode 1) or query (cov-mode 2) or both (cov-mode 0)
params.min_seq_id = 0.9 //[0.0-1.0] as the # of identical aligned residues divided by the number of aligned columns including internal gap columns: only matches above sequence identity __%
params.mmseqs_matrix = "${projectDir}/contrib/VTML10.out"
params.mmseqs_additional = "--clust-hash 1 --kmer-per-seq 100 --max-seqs 10000 -s 7.5 --alignment-mode 3 --sub-mat 'aa:${params.mmseqs_matrix},nucl:nucleotide.out' --seed-sub-mat 'aa:${params.mmseqs_matrix},nucl:nucleotide.out'" //additional mmseqs parameters (e.g. for more accurate recsults or higher speed)
//Note: MMSeqs2 automatically picks the optimal clustering strategy based on the coverage mode (--cov-mode 0 = set cover, --cov-mode 1,2 = greedy incremental); createdb-mode 1: softlink fasta

//fasta files parameters
params.fasta_prefix = "proteome_"
params.fasta_extension = ".fa"

//filtering parameters (minimum number of proteomes required to be present in cluster)
params.proteomes_threshold = 2 //default 2: keep all non singleton clusters

//config and constants
def set_config() {
    return [
    mmseqs_command: "easy-cluster", //alternative: easy-linclust
    label_clusters: "${projectDir}/scripts/label_clusters.py --batchsize 4000",
    filter_clusters: "${projectDir}/scripts/filter_clusters.py",
    extract_clusters: "${projectDir}/scripts/extract_clusters.py",
    align_clusters: "${projectDir}/scripts/dedup_align_fastafiles.py",
    presence_matrix: "${projectDir}/scripts/presence_matrix.py",
    fasta_indexer: "indexer.py -v -e '>' -i '^(.+\\|.*)\$' -r", //from ffdb.py
    fasta_extractor: "extractor.py -v", //from ffdb.py
    ]
}

process gather_fasta {
    // gather proteome fasta files present in given path and concatenate them
    maxForks params.maxForks.toInteger()
    cpus 1
    memory '1 GB'

    tag { name }

    input:
    tuple val(c), val(name), path(dir)

    output:
    tuple val(c), val(name), path(dir), path("${name}.fa"), optional: true

    script:
    """
    set -u
    set -e
    echo " .-- BEGUN ${task.process} for ${name} \$(date)"

    proteomefiles=\$(find ${dir}/ -type f -name "*.fa")
    proteomescount=\$(echo \${proteomefiles} | wc -w)
    if [ \$proteomescount -gt 0 ]; then
        echo " |-- Concatenating \${proteomescount} proteome fasta files... \$(date)"
        echo "\${proteomefiles}" #debug
        echo -n '' >${name}.fa
        for fastafile in \${proteomefiles}; do
            cat \${fastafile} >>${name}.fa
        done
        if [ ! -s "${name}.fa" ]; then
            rm ${name}.fa
            echo " *** ERROR, no fasta content found for ${name}"
        fi
    else
        echo " *** ERROR, no fasta files found under ${name}"
    fi

    echo " '-- ENDED ${task.process} \$(date)"
    """
}

process index_fasta {
    // index fasta file; note: requires ffdb.py
    maxForks params.maxForks.toInteger()
    cpus 1
    memory '1 GB'

    tag { name }

    input:
    tuple val(c), val(name), path(dir), path(fastaFile)

    output:
    tuple val(c), val(name), path(fastaFile), path("${fastaFile}.idx")

    script:
    """
    set -u
    set -e
    echo " .-- BEGUN ${task.process} for ${name} \$(date)"

    ${c.fasta_indexer} -f ${fastaFile} >${fastaFile}.idx

    echo " '-- ENDED ${task.process} \$(date)"
    """
}

process cluster_proteomes {
    // cluster proteome fasta files present in given path
    maxForks params.maxForks.toInteger()
    cpus params.clusterthreads
    memory '5 GB' //50

    tag { name }

    input:
    tuple val(c), val(name), path(dir), path(fastaFile)

    output:
    tuple val(c), val(name), path(dir), path("${name}_protein_cluster.tsv"), optional: true

    script:
    """
    set -u
    set -e
    echo " .-- BEGUN ${task.process} for ${name} \$(date)"

    echo " |-- Starting clustering job... \$(date)"
    command="${params.mmseqs_exe} ${c.mmseqs_command} ${fastaFile} ./_protein ./tmp --threads ${params.clusterthreads} --cov-mode ${params.covmode} --seq-id-mode ${params.seq_id_mode} --min-seq-id ${params.min_seq_id} -c ${params.coverage} --createdb-mode 1 ${params.mmseqs_additional}"

    echo " |-- \$command"
    eval \$command
    mv _protein_cluster.tsv ${name}_protein_cluster.tsv

    #cleanup
    rm -rf tmp/ _protein_all_seqs.fasta _protein_rep_seq.fasta
    echo " '-- ENDED ${task.process} \$(date)"
    """
}

process label_clusters {
    maxForks params.maxForks.toInteger()
    cpus params.labelthreads
    memory '3 GB' //30

    tag { name }

    input:
    tuple val(c), val(name), path(dir), path(clusterFile)

    output:
    tuple val(c), val(name), path("Labelled_protein_cluster.tsv")

    script:
    """
    set -u
    set -e
    echo " .-- BEGUN ${task.process} for ${name} \$(date)"

    ${c.label_clusters} --fasta_dir ${dir} --input_file ${clusterFile} --out_file Labelled_protein_cluster.tsv --prefix ${params.fasta_prefix} --extension ${params.fasta_extension} -t ${params.labelthreads} --sortlabels --uniq

    echo " '-- ENDED ${task.process} \$(date)"
    """
}

process filter_clusters {
    maxForks params.maxForks.toInteger()
    cpus 1
    memory '2 GB'

    tag { name }

    publishDir "${params.outdir}", mode: 'symlink'

    input:
    tuple val(c), val(name), path(labelledClustersFile)

    output:
    tuple val(c), val(name), path("${name}_clusters_m${params.proteomes_threshold}.tsv")

    script:
    """
    set -u
    set -e
    echo " .-- BEGUN ${task.process} for ${name} \$(date)"

    if [ "${params.proteomes_threshold}" -gt 1 ]; then
       echo " |-- keeping protein clusters with at least ${params.proteomes_threshold} proteomes"
    else
       echo " |-- keeping all protein clusters"
    fi
    ${c.filter_clusters} --input_file ${labelledClustersFile} --out_file ${name}_clusters_m${params.proteomes_threshold}.tsv --minproteomes ${params.proteomes_threshold} --strict

    echo " '-- ENDED ${task.process} \$(date)"
    """
}

process create_representatives {
    maxForks params.maxForks.toInteger()
    cpus params.extractthreads
    memory '2 GB' //10

    tag { name }

    publishDir "${params.outdir}", mode: 'symlink'

    input:
    tuple val(c), val(name), path(filteredClustersFile), path(fastaFile), path(indexFile)

    output:
    tuple val(c), val(name), path("${name}_repr.fa")

    script:
    """
    set -u
    set -e
    echo " .-- BEGUN ${task.process} for ${name} \$(date)"

    mkdir -p ${name}_clusters/
    tail -n +2 ${filteredClustersFile} | cut -f 5 >repr.ids
    ${c.fasta_extractor} -f ${fastaFile} -i ${indexFile} -l repr.ids >${name}_repr.fa
    rm repr.ids

    echo " '-- ENDED ${task.process} \$(date)"
    """
}

process extract_clusters {
    maxForks params.maxForks.toInteger()
    cpus params.extractthreads
    memory '2 GB' //10

    tag { name }

    publishDir "${params.outdir}", mode: 'symlink'

    input:
    tuple val(c), val(name), path(filteredClustersFile), path(fastaFile), path(indexFile)

    output:
    tuple val(c), val(name), path("${name}_clusters")

    script:
    """
    set -u
    set -e
    echo " .-- BEGUN ${task.process} for ${name} \$(date)"

    mkdir -p ${name}_clusters/
    ${c.extract_clusters} --input_file ${filteredClustersFile} --out_dir ${name}_clusters --single_file ${fastaFile} -t ${params.extractthreads} 

    echo " '-- ENDED ${task.process} \$(date)"
    """
}

process align_clusters {
    maxForks params.maxForks.toInteger()
    cpus params.alignthreads
    memory '2 GB' //400

    tag { name }

    publishDir "${params.outdir}", mode: 'symlink'

    input:
    tuple val(c), val(name), path(clustersDir)

    output:
    tuple val(c), val(name), path("${name}_alignments"), path("${name}_alnwidth.tsv")

    script:
    """
    set -u
    set -e
    echo " .-- BEGUN ${task.process} for ${name} \$(date)"

    mkdir -p ${name}_alignments/
    ${c.align_clusters} --fasta_dir ${clustersDir} --out_dir ${name}_alignments -t ${params.alignthreads} -at ${params.clustalothreads} -e .fa
    for i in ${name}_alignments/*.aln; do
        clusterid=\$(basename \$i);
        echo -en "\${clusterid%%.aln}\t"; awk '{if(NR==1) {print \$0} else {if(\$0 ~ /^>/) {print "\\n"\$0} else {printf \$0}}}' \$i | head -2 | wc -L
    done | sort -n -k1,1 >${name}_alnwidth.tsv


    echo " '-- ENDED ${task.process} \$(date)"
    """
}

process create_matrix {
    maxForks params.maxForks.toInteger()
    cpus 1
    memory '2 GB'

    tag { name }

    publishDir "${params.outdir}", mode: 'symlink'

    input:
    tuple val(c), val(name), path(filteredClustersFile)

    output:
    tuple val(c), val(name), path("${name}_m${params.proteomes_threshold}_matrix.tsv")

    script:
    """
    set -u
    set -e
    echo " .-- BEGUN ${task.process} for ${name} \$(date)"

    ${c.presence_matrix} -i ${filteredClustersFile} -p ${params.proteomesdir}/${name}.tsv -c -t -o ${name}_m${params.proteomes_threshold}_matrix.tsv

    echo " '-- ENDED ${task.process} \$(date)"
    """
}

workflow {
    config = set_config()
    println "NF:Cluster and Align (author: Giuseppe Insana)"
    println "-----"
    println "CONFIGURATION:"
    println "proteomesdir: ${params.proteomesdir}"
    println "outdir: ${params.outdir}"
    println "-----"
    println "covmode: ${params.covmode}"
    println "coverage: ${params.coverage}"
    println "seq_id_mode: ${params.seq_id_mode}"
    println "min_seq_id: ${params.min_seq_id}"
    println "-----"
    println "fasta_prefix: ${params.fasta_prefix}"
    println "fasta_extension: ${params.fasta_extension}"
    println "-----"
    println "proteomes_threshold: ${params.proteomes_threshold}"
    println "-----"

    def proteomesDir = file(params.proteomesdir)
    if( ! proteomesDir.exists() ) {
        log.error "Directory does not exist: ${proteomesDir}"
        System.exit(2)
    }
    def outDir = file(params.outdir)
    if( ! outDir.exists() ) {
        log.error "Directory does not exist: ${outDir}"
        System.exit(2)
    }

    //list steps that will be carried out
    if(params.align) {
        log.info "MAIN: cluster -> label -> filter -> extract -> align"
        log.info "PLUS: index, create_matrix"
    } else if(params.extract) {
        log.info "MAIN: cluster -> label -> filter -> extract"
        log.info "PLUS: index, create_matrix"
    } else {
        log.info "MAIN: cluster -> label -> filter"
        log.info "PLUS: create_matrix"
    }
    if(params.repr) {
        log.info "ALSO: create_repr"
    }
    println "-----"

    //1) channel with subfolders of the proteomesdir (e.g. proteomesdir/ecoli, proteomesdir/spneumo)
    channel.fromList( file(params.proteomesdir).listFiles()
           .findAll{ it.isDirectory() }
           .collect{ [it.name, file(it)] } ).set{proteome_dirs_ch}
    //proteome_dirs_ch.view() //debug

    proteome_dirs_ch = proteome_dirs_ch.map{ name, dir -> tuple(config, name, dir) }

    //2) concatenate proteomes' fasta files
    proteome_dirs_ch | gather_fasta

    //3) index input fasta files
    if(params.extract || params.align || params.repr) {
        gather_fasta.out | index_fasta
    }

    //4) cluster proteomes and label the clusters
    gather_fasta.out | cluster_proteomes | label_clusters

    //5) filter clusters to remove clusters with less proteomes than proteomes_threshold
    label_clusters.out | filter_clusters

    //6) create presence/absence/count matrix
    filter_clusters.out | create_matrix

    if(params.extract || params.align || params.repr) {
        def clusters_by_name = filter_clusters.out.map {
            c, name, labelledClustersFile ->
            tuple(name, labelledClustersFile)
        }
        def index_by_name = index_fasta.out.map {
            c, name, fastaFile, indexFile ->
            tuple(name, fastaFile, indexFile)
        }
        def extract_ch = clusters_by_name
            .combine(index_by_name, by: 0)
        extract_ch = extract_ch.map { name, labelledClustersFile, fastaFile, indexFile ->
            tuple(config, name, labelledClustersFile, fastaFile, indexFile)
        }

        //7) create fasta file with sequences of mmseqs cluster representatives
        if(params.repr) {
            extract_ch | create_representatives
        }

        //8) extract clusters
        if (params.extract || params.align) {
            extract_ch | extract_clusters
        }

        //9) align clusters
        if(params.align) {
            extract_clusters.out | align_clusters
        }
    }
}
