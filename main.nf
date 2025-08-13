#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

//modules

process FASTQC {
    tag "$meta.id"
    label 'process_low'
    
    conda "bioconda::fastqc=0.12.1"
    container "quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"
    
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip") , emit: zip
    path "versions.yml"           , emit: versions

    script:
    """
    fastqc --quiet --threads $task.cpus $reads

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$(fastqc --version | sed -e "s/FastQC v//g")
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_fastqc.html
    touch ${meta.id}_fastqc.zip

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$(fastqc --version | sed -e "s/FastQC v//g")
    END_VERSIONS
    """
}

process STAR_ALIGN {
    tag "$meta.id"
    label 'process_high'
    
    conda "bioconda::star=2.7.11a bioconda::samtools=1.19.2"
    container "quay.io/biocontainers/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:1df389393721fc66f3fd8778ad938ac711951107-0"
    
    publishDir "${params.outdir}/star", mode: 'copy'

    input:
    tuple val(meta), path(reads)
    path index
    path gtf

    output:
    tuple val(meta), path('*.bam')        , emit: bam
    tuple val(meta), path('*.out')        , emit: log
    tuple val(meta), path('*SJ.out.tab')  , emit: sj
    path "versions.yml"                   , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    STAR \\
        --genomeDir $index \\
        --readFilesIn $reads \\
        --runThreadN $task.cpus \\
        --outFileNamePrefix ${prefix}. \\
        --outSAMtype BAM SortedByCoordinate \\
        --sjdbGTFfile $gtf \\
        --readFilesCommand zcat \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.Aligned.sortedByCoord.out.bam
    touch ${prefix}.Log.out
    touch ${prefix}.Log.final.out
    touch ${prefix}.SJ.out.tab

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}

process STAR_GENOMEGENERATE {
    tag "star_genomegenerate"
    label 'process_high'
    
    conda "bioconda::star=2.7.11a bioconda::samtools=1.19.2"
    container "quay.io/biocontainers/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:1df389393721fc66f3fd8778ad938ac711951107-0"
    
    publishDir "${params.outdir}/star_index", mode: 'copy'

    input:
    path fasta
    path gtf

    output:
    path "star_index"     , emit: index
    path "versions.yml"   , emit: versions

    script:
    def args = task.ext.args ?: ''
    def memory_gb = task.memory ? "--limitGenomeGenerateRAM ${task.memory.toBytes() - 100000000}" : ''
    """
    mkdir star_index

    STAR \\
        --runMode genomeGenerate \\
        --genomeDir star_index/ \\
        --genomeFastaFiles $fasta \\
        --sjdbGTFfile $gtf \\
        --runThreadN $task.cpus \\
        --sjdbOverhang 100 \\
        $memory_gb \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
    END_VERSIONS
    """

    stub:
    """
    mkdir star_index
    touch star_index/chrLength.txt
    touch star_index/chrNameLength.txt
    touch star_index/chrName.txt
    touch star_index/chrStart.txt
    touch star_index/Genome
    touch star_index/genomeParameters.txt
    touch star_index/SA
    touch star_index/SAindex

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
    END_VERSIONS
    """
}

process FEATURECOUNTS {
    tag "$meta.id"
    label 'process_medium'
    
    conda "bioconda::subread=2.0.6"
    container "quay.io/biocontainers/subread:2.0.6--he4a0461_0"
    
    publishDir "${params.outdir}/featurecounts", mode: 'copy'

    input:
    tuple val(meta), path(bam)
    path gtf

    output:
    tuple val(meta), path("*featureCounts.txt")        , emit: counts
    tuple val(meta), path("*featureCounts.txt.summary"), emit: summary
    path "versions.yml"                                , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    featureCounts \\
        $args \\
        -T $task.cpus \\
        -a $gtf \\
        -o ${prefix}.featureCounts.txt \\
        $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        subread: \$(featureCounts -v 2>&1 | grep -E -o 'v[0-9]+\\.[0-9]+\\.[0-9]+' | sed 's/v//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.featureCounts.txt
    touch ${prefix}.featureCounts.txt.summary

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        subread: \$(featureCounts -v 2>&1 | grep -E -o 'v[0-9]+\\.[0-9]+\\.[0-9]+' | sed 's/v//')
    END_VERSIONS
    """
}
// Parameters
params.input = 'samplesheet.csv'
params.fasta = null
params.gtf = null
params.star_index = null
params.outdir = 'results'

// Help message
def helpMessage() {
    log.info"""
    Usage:
    nextflow run main.nf --input samplesheet.csv --fasta genome.fa --gtf annotation.gtf

    Required arguments:
      --input         Path to samplesheet CSV file
      --gtf           Path to GTF annotation file

    Genome input (choose one):
      --fasta         Path to genome FASTA file (will build STAR index)
      --star_index    Path to pre-built STAR genome index directory

    Optional arguments:
      --outdir        Output directory (default: results)

    Examples:
      # Build index from FASTA
      nextflow run main.nf --input samplesheet.csv --fasta genome.fa --gtf annotation.gtf

      # Use pre-built index
      nextflow run main.nf --input samplesheet.csv --star_index /path/to/index --gtf annotation.gtf
    """.stripIndent()
}

// Show help message if requested
if (params.help) {
    helpMessage()
    exit 0
}

// Validate required parameters
if (!params.gtf) {
    error "Please provide a GTF file with --gtf"
}

if (!params.fasta && !params.star_index) {
    error "Please provide either a genome FASTA file (--fasta) or a pre-built STAR index (--star_index)"
}

if (params.fasta && params.star_index) {
    error "Please provide either --fasta OR --star_index, not both"
}

workflow {
    // Read samplesheet and create input channel
    ch_input = Channel
        .fromPath(params.input)
        .splitCsv(header: true)
        .map { row ->
            def meta = [id: row.sample, single_end: true]
            [meta, file(row.fastq_1)]
        }

    // Run FastQC
    FASTQC(ch_input)

    // Handle STAR index - either build it or use existing
    if (params.fasta) {
        // Build STAR index from FASTA
        STAR_GENOMEGENERATE(
            file(params.fasta),
            file(params.gtf)
        )
        ch_star_index = STAR_GENOMEGENERATE.out.index
    } else {
        // Use pre-built index
        ch_star_index = Channel.value(file(params.star_index))
    }

    // Run STAR alignment
    STAR_ALIGN(
        ch_input,
        ch_star_index,
        file(params.gtf)
    )

    // Run featureCounts
    FEATURECOUNTS(
        STAR_ALIGN.out.bam,
        file(params.gtf)
    )
}