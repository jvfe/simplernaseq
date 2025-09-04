#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

//modules
process GUNZIP {
    tag "$archive"
    label 'process_single'

    conda "conda-forge::sed=4.7"
    container "quay.io/nf-core/ubuntu:20.04"

    input:
    tuple val(meta), path(archive)

    output:
    tuple val(meta), path("$gunzip"), emit: gunzip
    path "versions.yml"             , emit: versions

    script:
    def args = task.ext.args ?: ''
    gunzip = archive.toString() - '.gz'
    """
    # Not calling gunzip itself because it can be very slow
    # for large files. Using bash and zcat instead.
    zcat $args $archive > $gunzip
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gunzip: \$(echo \$(gunzip --version 2>&1) | sed 's/^.*(gzip) //; s/ Copyright.*\$//')
    END_VERSIONS
    """

    stub:
    gunzip = archive.toString() - '.gz'
    """
    touch $gunzip
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gunzip: \$(echo \$(gunzip --version 2>&1) | sed 's/^.*(gzip) //; s/ Copyright.*\$//')
    END_VERSIONS
    """
}

process GUNZIP_GTF {
    tag "$archive"
    label 'process_single'

    conda "conda-forge::sed=4.7"
    container "quay.io/nf-core/ubuntu:20.04"

    input:
    tuple val(meta), path(archive)

    output:
    tuple val(meta), path("$gunzip"), emit: gunzip
    path "versions.yml"             , emit: versions

    script:
    def args = task.ext.args ?: ''
    gunzip = archive.toString() - '.gz'
    """
    # Not calling gunzip itself because it can be very slow
    # for large files. Using bash and zcat instead.
    zcat $args $archive > $gunzip

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gunzip: \$(echo \$(gunzip --version 2>&1) | sed 's/^.*(gzip) //; s/ Copyright.*\$//')
    END_VERSIONS
    """

    stub:
    gunzip = archive.toString() - '.gz'
    """
    touch $gunzip
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gunzip: \$(echo \$(gunzip --version 2>&1) | sed 's/^.*(gzip) //; s/ Copyright.*\$//')
    END_VERSIONS
    """
}


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
    
    conda "bioconda::hisat2=2.2.1 bioconda::samtools=1.16.1"
    container "quay.io/biocontainers/mulled-v2-a97e90b3b802d1da3d6958e0867610c718cb5eb1:2cdf6bf1e92acbeb9b2834b1c58754167173a410-0"
    
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
    
    conda "bioconda::hisat2=2.2.1 bioconda::samtools=1.16.1"
    container "quay.io/biocontainers/mulled-v2-a97e90b3b802d1da3d6958e0867610c718cb5eb1:2cdf6bf1e92acbeb9b2834b1c58754167173a410-0"
    
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

process CUSTOM_DUMPSOFTWAREVERSIONS {
    label 'process_single'

    conda "bioconda::multiqc=1.14"
    container "quay.io/biocontainers/multiqc:1.14--pyhdfd78af_0"

    input:
    path versions

    output:
    path "software_versions.yml"       , emit: yml
    path "software_versions_mqc.yml"   , emit: mqc_yml

    script:
    """
    #!/usr/bin/env python3

    import yaml
    import platform
    from textwrap import dedent

    def _make_versions_html(versions):
        html = [
            dedent(
                '''\\
                <style>
                #nf-core-versions tbody:nth-child(even) {
                    background-color: #f2f2f2;
                }
                </style>
                <table class="table" style="width:100%" id="nf-core-versions">
                    <thead>
                        <tr>
                            <th> Process Name </th>
                            <th> Software </th>
                            <th> Version  </th>
                        </tr>
                    </thead>
                '''
            )
        ]
        for process, tmp_versions in sorted(versions.items()):
            html.append("<tbody>")
            for i, (tool, version) in enumerate(sorted(tmp_versions.items())):
                html.append(
                    dedent(
                        f'''\\
                        <tr>
                            <td><samp>{process if (i == 0) else ""}</samp></td>
                            <td><samp>{tool}</samp></td>
                            <td><samp>{version}</samp></td>
                        </tr>
                        '''
                    )
                )
            html.append("</tbody>")
        html.append("</table>")
        return "".join(html)

    with open("$versions") as f:
        versions_this_module = yaml.safe_load(f)

    # Aggregate versions by the module name:
    versions_by_module = {}
    for process_name, process_versions in versions_this_module.items():
        module_name = process_name.split(":")[-1]
        try:
            if versions_by_module[module_name] != process_versions:
                raise AssertionError(
                    "We assume that software versions are the same between all modules. "
                    "If you see this error there are probably subworkflows that include "
                    "the same module with different versions."
                )
        except KeyError:
            versions_by_module[module_name] = process_versions

    versions_by_module["Workflow"] = {
        "Nextflow": "$workflow.nextflow.version",
        "$workflow.manifest.name": "$workflow.manifest.version"
    }

    versions_mqc = {
        'id': 'software_versions',
        'section_name': 'Software Versions',
        'section_href': 'https://github.com/nf-core/tools',
        'plot_type': 'html',
        'description': 'are collected at run time from the software output.',
        'data': _make_versions_html(versions_by_module)
    }

    with open("software_versions.yml", 'w') as f:
        yaml.dump(versions_by_module, f, default_flow_style=False)
    with open("software_versions_mqc.yml", 'w') as f:
        yaml.dump(versions_mqc, f, default_flow_style=False)

    with open('versions.yml', 'w') as f:
        yaml.dump(versions_by_module, f, default_flow_style=False)
    """

    stub:
    """
    touch software_versions.yml
    touch software_versions_mqc.yml
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

    // Run FastQC
    FASTQC(ch_input)

    // Handle GTF file - decompress if needed
    ch_gtf = Channel.fromPath(params.gtf)
        .map { gtf_file ->
            def meta = [id: 'gtf']
            [meta, gtf_file]
        }
    
    if (params.gtf.endsWith('.gz')) {
        GUNZIP_GTF(ch_gtf)
        ch_gtf_final = GUNZIP_GTF.out.gunzip.map { meta, file -> file }
        ch_gunzip_gtf_versions = GUNZIP_GTF.out.versions
    } else {
        ch_gtf_final = ch_gtf.map { meta, file -> file }
        ch_gunzip_gtf_versions = Channel.empty()
    }

    // Handle HISAT2 index - either build it or use existing
    if (params.fasta) {
        // Handle FASTA file - decompress if needed
        ch_fasta = Channel.fromPath(params.fasta)
            .map { fasta_file ->
                def meta = [id: 'fasta']
                [meta, fasta_file]
            }
        
        if (params.fasta.endsWith('.gz')) {
            GUNZIP(ch_fasta)
            ch_fasta_final = GUNZIP.out.gunzip.map { meta, file -> file }
            ch_gunzip_fasta_versions = GUNZIP.out.versions
        } else {
            ch_fasta_final = ch_fasta.map { meta, file -> file }
            ch_gunzip_fasta_versions = Channel.empty()
        }

        // Build HISAT2 index from FASTA
        HISAT2_BUILD(
            ch_fasta_final,
            ch_gtf_final,
            []
        )
        ch_hisat2_index = HISAT2_BUILD.out.index.first()
        ch_hisat2_build_versions = HISAT2_BUILD.out.versions
    } else {
        // Use pre-built index
        ch_hisat2_index = Channel.value(file(params.hisat2_index))
        ch_hisat2_build_versions = Channel.empty()
        ch_gunzip_fasta_versions = Channel.empty()
    }

    // Run HISAT2 alignment
    HISAT2_ALIGN(
        ch_input,
        ch_hisat2_index,
        ch_gtf_final.first()
    )

    // Run featureCounts
    FEATURECOUNTS(
        HISAT2_ALIGN.out.bam,
        ch_gtf_final.first()
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
