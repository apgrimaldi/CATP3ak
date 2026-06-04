process PICARD_MARKDUPLICATES {
    tag "$meta.id"
    label 'process_medium'
    container 'quay.io/biocontainers/picard:2.27.4--hdfd78af_0'

    input:
    tuple val(meta), path(bam)
    path fasta  // File singolo del genoma
    path fai    // File singolo dell'indice .fai

    output:
    tuple val(meta), path("*.removed.bam")       , emit: bam
    tuple val(meta), path("*.{bai,bam.bai}")     , emit: bai, optional: true
    tuple val(meta), path("*.metrics.txt")       , emit: metrics
    path  "versions.yml"                         , emit: versions

    script:
    def args = task.ext.args ?: '--REMOVE_DUPLICATES true --ASSUME_SORTED true --VALIDATION_STRINGENCY LENIENT --CREATE_INDEX true'
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    
    def reference = fasta ? "--REFERENCE_SEQUENCE ${fasta}" : ""
    
    
    def avail_mem = 3072
    if (task.memory) {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }

    """
    picard \\
        -Xmx${avail_mem}M \\
        MarkDuplicates \\
        $args \\
        --INPUT $bam \\
        --OUTPUT ${prefix}.removed.bam \\
        --METRICS_FILE ${prefix}.metrics.txt \\
        $reference

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard MarkDuplicates --version 2>&1 | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """
}
