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
    path "versions.yml"                     , emit: versions

    script:
    def args = task.ext.args ?: ''
    
    """
    # Esegue genomeDistribution in modalità dettagliata (-d)
    genomeDistribution -d -i ${peaks} -g ${gtf} -o ${meta.id}_distribution ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqcode: \$(echo "base")
    END_VERSIONS
    """
}
