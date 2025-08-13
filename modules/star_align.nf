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