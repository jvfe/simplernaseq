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