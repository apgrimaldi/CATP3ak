process DEEPTOOLS_PLOTHEATMAP {
    tag "$meta.id"
    label 'process_medium'

    container 'biocontainers/deeptools:3.5.5--pyhdfd78af_0'

    input:
    tuple val(meta), path(matrix)

    output:
    tuple val(meta), path("*.png"), emit: png
    tuple val(meta), path("*.pdf"), emit: pdf
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """

    plotHeatmap \\
        --matrixFile ${matrix} \\
        --outFileName ${prefix}.heatmap.png \\
        $args

    
    plotHeatmap \\
        --matrixFile ${matrix} \\
        --outFileName ${prefix}.heatmap.pdf \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: \$(plotHeatmap --version | sed -e "s/plotHeatmap //g")
    END_VERSIONS
    """
}
