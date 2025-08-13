#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Include modules
include { FASTQC } from './modules/fastqc.nf'
include { STAR_GENOMEGENERATE } from './modules/star_genomegenerate.nf'  // Fixed typo too
include { STAR_ALIGN } from './modules/star_align.nf'
include { FEATURECOUNTS } from './modules/featurecounts.nf'

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