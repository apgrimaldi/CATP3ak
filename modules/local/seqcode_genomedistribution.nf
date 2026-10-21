process SEQCODE_GENOMEDISTRIBUTION {
    tag "$meta.id"
    label 'process_medium'
    
    container 'eblancocrg/seqcode:base'

    publishDir "${params.outdir}/14_seqcode_distribution", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(peaks)
    path gtf

    output:
    tuple val(meta), path("*_distribution*"), emit: results

    script:
    def args = task.ext.args ?: ''
    
    """
    # 1. Cerca il percorso esatto dell'eseguibile dentro il container
    CMD_PATH=\$(find / -name "genomeDistribution" -type f -executable 2>/dev/null | grep "bin" | head -n 1)
    
    # Fallback se find non dovesse trovarlo (percorso tipico di installazione)
    if [ -z "\$CMD_PATH" ]; then
        CMD_PATH="/seqcode/bin/genomeDistribution"
    fi

    # 2. Esegue il tool usando il percorso assoluto e i parametri
    \$CMD_PATH -d -i ${peaks} -g ${gtf} -o ${meta.id}_distribution ${args}
    """
}
