process MACS3_CHIP_NARROW {
    tag "$meta.id"
    label 'process_medium'
    container 'quay.io/biocontainers/macs3:3.0.1--py311h0152c62_3'

    input:
    tuple val(meta), path(ip_bam), path(control_bam)
    val gsize

    output:
    tuple val(meta), path("*.narrowPeak")       , emit: peaks
    tuple val(meta), path("*.xls")              , emit: xls
    tuple val(meta), path("*.narrow_counts.txt"), emit: count_narrow
    path "versions.yml"                         , emit: versions

    script:
    def prefix   = "${meta.id}_narrow"
    def format   = meta.single_end ? 'BAM' : 'BAMPE'
    def args_control = control_bam ? "-c $control_bam" : ""

    """
    macs3 callpeak \\
        -t $ip_bam \\
        $args_control \\
        -f $format \\
        -g $gsize \\
        -n $prefix \\
        --qvalue 0.05

    
    PEAK_FILE=\$(ls *.narrowPeak 2>/dev/null | head -n 1)
    if [ -n "\$PEAK_FILE" ] && [ -f "\$PEAK_FILE" ]; then
        count=\$(wc -l < "\$PEAK_FILE")
    else
        count=0
    fi

    echo -e "Sample\\tNarrow_Peaks" > ${prefix}.narrow_counts.txt
    echo -e "${meta.id}\\t\$count" >> ${prefix}.narrow_counts.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        macs3: \$(macs3 --version | sed 's/macs3 //g')
    END_VERSIONS
    """
}
