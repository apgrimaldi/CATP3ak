process BOWTIE2_BUILD {
    tag "$fasta"
    label 'process_high'
    container 'quay.io/biocontainers/bowtie2:2.5.5--ha27dd3b_0' 

    input:
    path fasta

    output:
    path "*.bt2*", emit: index
    path "versions.yml", emit: versions

    script:
    """
    bowtie2-build --threads $task.cpus $fasta genome_index
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(echo \$(bowtie2 --version 2>&1) | sed 's/^.*bowtie2-align-s version //; s/ .*\$//')
    END_VERSIONS
    """
}
