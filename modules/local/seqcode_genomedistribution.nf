process SEQCODE_GENOMEDISTRIBUTION {
    tag "$meta.id"
    label 'process_medium'

    container 'eblancocrg/seqcode:base'

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
    
    CMD_PATH=\$(find / -name "genomeDistribution" -type f -executable 2>/dev/null | grep "bin" | head -n 1)
    if [ -z "\$CMD_PATH" ]; then
        CMD_PATH="/seqcode/bin/genomeDistribution"
    fi

    
    REFGENE_UNCOMPRESSED="refGene.txt"
    if [[ "${refgene}" == *.gz ]]; then
        gunzip -c "${refgene}" > \$REFGENE_UNCOMPRESSED
    else
        cp "${refgene}" \$REFGENE_UNCOMPRESSED
    fi

    
    # <chrom_info> <refGene.txt> <peaks.bed> <name>
    \$CMD_PATH \\
        $args \\
        ${chrom_sizes} \\
        \$REFGENE_UNCOMPRESSED \\
        ${peak} \\
        ${meta.id}_seqcode

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqcode: "1.0"
    END_VERSIONS
    """
}
