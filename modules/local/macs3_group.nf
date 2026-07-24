process MACS3_GROUP {
    tag "$meta.id"
    label 'process_high'
    container 'quay.io/biocontainers/macs3:3.0.1--py311h0152c62_3'

    input:
    tuple val(meta), path(ip_bams), path(control_bams)
    val macs_gsize

    output:
    tuple val(meta), path("*_peaks.narrowPeak"), emit: peaks
    tuple val(meta), path("*_peaks.xls")       , emit: xls
    path "*_peaks.xls"                         , emit: count_narrow
    path "versions.yml"                        , emit: versions

    script:
    def args = task.ext.args ?: '--nomodel --shift -100 --extsize 200 -B'

    """
    macs3 callpeak \\
        -t ${ip_bams} \\
        -c ${control_bams} \\
        -f BAM \\
        -g ${macs_gsize} \\
        -n ${meta.id} ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        macs3: \$(macs3 --version | sed -e "s/macs3 //g")
    END_VERSIONS
    """
}
