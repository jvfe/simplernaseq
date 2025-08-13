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

process HISAT2_ALIGN {
    tag "$meta.id"
    label 'process_high'
    
    conda "bioconda::hisat2=2.2.1 bioconda::samtools=1.19.2"
    container "quay.io/biocontainers/mulled-v2-a97e90b3b802d1da3d6958e0867610c718cb5eb1:2880dd9d8ad0a7b221d4ebc6e532a6ac915b6e90-0"
    
    publishDir "${params.outdir}/hisat2", mode: 'copy'

    input:
    tuple val(meta), path(reads)
    path index
    path gtf

    output:
    tuple val(meta), path('*.bam')        , emit: bam
    tuple val(meta), path('*.log')        , emit: log
    path "versions.yml"                   , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def strandedness = ''
    if (meta.strandedness == 'forward') {
        strandedness = meta.single_end ? '--rna-strandness F' : '--rna-strandness FR'
    } else if (meta.strandedness == 'reverse') {
        strandedness = meta.single_end ? '--rna-strandness R' : '--rna-strandness RF'
    }
    def seq_center = meta.seq_center ? "--rg-id ${prefix} --rg SM:${prefix} --rg CN:${meta.seq_center}" : "--rg-id ${prefix} --rg SM:${prefix}"
    """
    INDEX=`find -L ./ -name "*.1.ht2" | sed 's/\\.1\\.ht2\$//'`
    
    hisat2 \\
        -x \$INDEX \\
        $strandedness \\
        $seq_center \\
        --threads $task.cpus \\
        $args \\
        -U $reads \\
        --summary-file ${prefix}.hisat2.summary.log \\
        | samtools sort --threads $task.cpus -o ${prefix}.bam -

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hisat2: \$(echo \$(hisat2 --version 2>&1) | sed 's/^.*hisat2-align-s version //; s/ .*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    touch ${prefix}.hisat2.summary.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hisat2: \$(echo \$(hisat2 --version 2>&1) | sed 's/^.*hisat2-align-s version //; s/ .*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}

process HISAT2_BUILD {
    tag "hisat2_build"
    label 'process_high'
    
    conda "bioconda::hisat2=2.2.1 bioconda::samtools=1.19.2"
    container "quay.io/biocontainers/mulled-v2-a97e90b3b802d1da3d6958e0867610c718cb5eb1:2880dd9d8ad0a7b221d4ebc6e532a6ac915b6e90-0"
    
    publishDir "${params.outdir}/hisat2_index", mode: 'copy'

    input:
    path fasta
    path gtf
    path splicesites

    output:
    path "hisat2_index"     , emit: index
    path "versions.yml"     , emit: versions

    script:
    def args = task.ext.args ?: ''
    def avail_mem = 3
    if (!task.memory) {
        log.info '[HISAT2 index build] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    def ss = splicesites ? "--ss $splicesites" : ""
    """
    mkdir hisat2_index
    hisat2-build \\
        -p $task.cpus \\
        $ss \\
        $args \\
        $fasta \\
        hisat2_index/${fasta.baseName}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hisat2: \$(echo \$(hisat2 --version 2>&1) | sed 's/^.*hisat2-align-s version //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    """
    mkdir hisat2_index
    touch hisat2_index/${fasta.baseName}.1.ht2
    touch hisat2_index/${fasta.baseName}.2.ht2
    touch hisat2_index/${fasta.baseName}.3.ht2
    touch hisat2_index/${fasta.baseName}.4.ht2
    touch hisat2_index/${fasta.baseName}.5.ht2
    touch hisat2_index/${fasta.baseName}.6.ht2
    touch hisat2_index/${fasta.baseName}.7.ht2
    touch hisat2_index/${fasta.baseName}.8.ht2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hisat2: \$(echo \$(hisat2 --version 2>&1) | sed 's/^.*hisat2-align-s version //; s/ .*\$//')
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

process MULTIQC {
    label 'process_single'
    
    conda "bioconda::multiqc=1.21"
    container "quay.io/biocontainers/multiqc:1.21--pyhdfd78af_0"
    
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path '*'
    val multiqc_title

    output:
    path "*multiqc_report.html", emit: report
    path "*_data"              , emit: data
    path "*_plots"             , optional:true, emit: plots
    path "versions.yml"        , emit: versions

    script:
    def args = task.ext.args ?: ''
    def custom_config = params.multiqc_config ? "--config $params.multiqc_config" : ''
    def extra_config = params.multiqc_methods_description ? "--cl-config 'custom_plot_config: {mqc_methods_description: \"${params.multiqc_methods_description}\"}'" : ''
    def logo = params.multiqc_logo ? "--cl-config 'custom_logo: \"${params.multiqc_logo}\"'" : ''
    """
    multiqc \\
        --force \\
        $args \\
        $custom_config \\
        $extra_config \\
        $logo \\
        --title "$multiqc_title" \\
        .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """

    stub:
    """
    touch ${multiqc_title}_multiqc_report.html
    mkdir ${multiqc_title}_data
    mkdir ${multiqc_title}_plots

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """
}

// Parameters
params.input = 'samplesheet.csv'
params.fasta = null
params.gtf = null
params.hisat2_index = null
params.outdir = 'results'
params.multiqc_config = null
params.multiqc_title = null
params.multiqc_logo = null
params.multiqc_methods_description = null

// Help message
def helpMessage() {
    log.info"""
    Usage:
    nextflow run main.nf --input samplesheet.csv --fasta genome.fa --gtf annotation.gtf

    Required arguments:
      --input         Path to samplesheet CSV file
      --gtf           Path to GTF annotation file

    Genome input (choose one):
      --fasta         Path to genome FASTA file (will build HISAT2 index)
      --hisat2_index  Path to pre-built HISAT2 genome index directory

    Optional arguments:
      --outdir        Output directory (default: results)
      --multiqc_title Custom title for MultiQC report
      --multiqc_config Path to MultiQC config file
      --multiqc_logo  Path to logo for MultiQC report

    Examples:
      # Build index from FASTA
      nextflow run main.nf --input samplesheet.csv --fasta genome.fa --gtf annotation.gtf

      # Use pre-built index
      nextflow run main.nf --input samplesheet.csv --hisat2_index /path/to/index --gtf annotation.gtf
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

if (!params.fasta && !params.hisat2_index) {
    error "Please provide either a genome FASTA file (--fasta) or a pre-built HISAT2 index (--hisat2_index)"
}

if (params.fasta && params.hisat2_index) {
    error "Please provide either --fasta OR --hisat2_index, not both"
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

    // Initialize versions channel
    ch_versions = Channel.empty()

    // Run FastQC
    FASTQC(ch_input)
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    // Handle HISAT2 index - either build it or use existing
    if (params.fasta) {
        // Build HISAT2 index from FASTA
        HISAT2_BUILD(
            file(params.fasta),
            file(params.gtf),
            []
        )
        ch_hisat2_index = HISAT2_BUILD.out.index
        ch_versions = ch_versions.mix(HISAT2_BUILD.out.versions)
    } else {
        // Use pre-built index
        ch_hisat2_index = Channel.value(file(params.hisat2_index))
    }

    // Run HISAT2 alignment
    HISAT2_ALIGN(
        ch_input,
        ch_hisat2_index,
        file(params.gtf)
    )
    ch_versions = ch_versions.mix(HISAT2_ALIGN.out.versions.first())

    // Run featureCounts
    FEATURECOUNTS(
        HISAT2_ALIGN.out.bam,
        file(params.gtf)
    )
    ch_versions = ch_versions.mix(FEATURECOUNTS.out.versions.first())

    // Collect all outputs for MultiQC
    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.html.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(HISAT2_ALIGN.out.log.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(FEATURECOUNTS.out.counts.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(FEATURECOUNTS.out.summary.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_versions.collect().ifEmpty([]))

    // Set MultiQC title
    multiqc_title = params.multiqc_title ?: 'RNA-seq Analysis Report'

    // Run MultiQC
    MULTIQC(
        ch_multiqc_files.collect(),
        multiqc_title
    )
}
