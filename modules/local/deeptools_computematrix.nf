process DEEPTOOLS_COMPUTEMATRIX {
    tag "$meta.id"
    label 'process_medium'

    container 'quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0'

    input:
    tuple val(meta), path(bigwig)
    path  bed

    output:
    tuple val(meta), path("*.mat.gz") , emit: matrix
    tuple val(meta), path("*.mat.tab"), emit: table
    path  "versions.yml"              , emit: versions

    script:
    def args   = task.ext.args ?: '--referencePoint TSS -b 3000 -a 3000 --binSize 50 --skipZeros'
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    computeMatrix reference-point \\
        $args \\
        --regionsFileName $bed \\
        --scoreFileName $bigwig \\
        --outFileName ${prefix}.computeMatrix.mat.gz \\
        --outFileNameMatrix ${prefix}.computeMatrix.vals.mat.tab \\
        --numberOfProcessors ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: \$(computeMatrix --version | sed -e "s/computeMatrix //g")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.computeMatrix.mat.gz
    touch ${prefix}.computeMatrix.vals.mat.tab

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: 3.5.5
    END_VERSIONS
    """
}
