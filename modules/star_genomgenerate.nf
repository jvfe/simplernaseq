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