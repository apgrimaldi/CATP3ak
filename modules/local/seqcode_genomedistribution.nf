process SEQCODE_GENOMEDISTRIBUTION {
    tag "$meta.id"
    label 'process_medium'

    container 'quay.io/biocontainers/seqcode:...' // (mantieni il container che già usi)

    input:
    tuple val(meta), path(peak)
    path refgene
    path chrom_sizes

    output:
    tuple val(meta), path("*.{txt,pdf,png}"), emit: results
    path "versions.yml"                   , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    
    REFGENE_UNCOMPRESSED="${refgene}"
    if [[ "${refgene}" == *.gz ]]; then
        gunzip -c "${refgene}" > refgene_unpacked.txt
        REFGENE_UNCOMPRESSED="refgene_unpacked.txt"
    fi

    
    seqcode genomedistribution \\
        $args \\
        --peaks ${peak} \\
        --refgene \$REFGENE_UNCOMPRESSED \\
        --chrom_sizes ${chrom_sizes} \\
        --output ${meta.id}_seqcode

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqcode: \$(seqcode --version 2>&1 | head -n 1 | sed 's/v//g')
    END_VERSIONS
    """
}
