process SEQCODE_GENOMEDISTRIBUTION {
    tag "$meta.id"
    label 'process_medium'
    
    container 'eblancocrg/seqcode:base'

    publishDir "${params.outdir}/14_seqcode_distribution", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(peaks)
    path annotation
    path chrom_sizes

    output:
    tuple val(meta), path("*"), emit: results

    script:
    def args = task.ext.args ?: ''
    
    """
    
    CMD_PATH=\$(find / -name "genomeDistribution" -type f -executable 2>/dev/null | grep "bin" | head -n 1)
    
    if [ -z "\$CMD_PATH" ]; then
        CMD_PATH="/seqcode/bin/genomeDistribution"
    fi

    
    \$CMD_PATH -d ${args} ${chrom_sizes} ${annotation} ${peaks} ${meta.id}_distribution
    """
}
