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