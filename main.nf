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

    // Handle HISAT2 index - either build it or use existing
    if (params.fasta) {
        // Build HISAT2 index from FASTA
        HISAT2_BUILD(
            file(params.fasta),
            file(params.gtf),
            []
        )
        ch_hisat2_index = HISAT2_BUILD.out.index
        ch_hisat2_build_versions = HISAT2_BUILD.out.versions
    } else {
        // Use pre-built index
        ch_hisat2_index = Channel.value(file(params.hisat2_index))
        ch_hisat2_build_versions = Channel.empty()
    }

    // Run HISAT2 alignment
    HISAT2_ALIGN(
        ch_input,
        ch_hisat2_index,
        file(params.gtf)
    )

    // Run featureCounts
    FEATURECOUNTS(
        HISAT2_ALIGN.out.bam,
        file(params.gtf)
    )

    // Collect versions from all processes
    ch_versions = Channel.empty()
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())
    ch_versions = ch_versions.mix(HISAT2_ALIGN.out.versions.first())
    ch_versions = ch_versions.mix(FEATURECOUNTS.out.versions.first())
    ch_versions = ch_versions.mix(ch_hisat2_build_versions)

    // Create consolidated versions file
    CUSTOM_DUMPSOFTWAREVERSIONS(
        ch_versions.collectFile(name: 'collated_versions.yml')
    )

    // Collect all outputs for MultiQC
    ch_multiqc_files = Channel.empty()
    
    // Add FastQC outputs
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    
    // Add HISAT2 alignment logs
    ch_multiqc_files = ch_multiqc_files.mix(HISAT2_ALIGN.out.log.collect{it[1]})
    
    // Add featureCounts outputs
    ch_multiqc_files = ch_multiqc_files.mix(FEATURECOUNTS.out.summary.collect{it[1]})
    
    // Add consolidated versions
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml)

    // Set MultiQC title
    multiqc_title = params.multiqc_title ?: 'RNA-seq Analysis Report'

    // Run MultiQC
    MULTIQC(
        ch_multiqc_files.collect(),
        multiqc_title
    )
}
