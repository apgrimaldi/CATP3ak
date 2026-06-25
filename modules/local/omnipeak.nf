process OMNIPEAK {
    tag "${meta.id}"
    label 'process_high'
    container 'biohaz/omnipeak:latest'

    input:
    tuple val(meta), path(ip_bam), path(control_bam)

    output:
    tuple val(meta), path("${meta.id}_peaks.bed"), emit: peaks
    path "*_omnipeak_mqc.txt"                    , emit: counts_mqc
    path "versions.yml"                          , emit: versions

    script:
    def has_control = control_bam.name != 'null' && control_bam.toString() != ''
    def control_arg = has_control ? "-c ${control_bam}" : ""
    
    """
    java --add-modules=jdk.incubator.vector -Xmx8G -jar /home/omnipeak/build/libs/omnipeak-1.5.build.jar analyze \\
        -t ${ip_bam} \\
        ${control_arg} \\
        --cs /home/omnipeak/chrom.sizes \\
        -p ${meta.id}_peaks.bed

    PEAK_COUNT=\$(wc -l < ${meta.id}_peaks.bed)
    cat <<EOF > ${meta.id}_omnipeak_mqc.txt
# id: 'omnipeak_counts'
# section_name: 'Omnipeak: Peaks Identified'
# plot_type: 'bargraph'
Sample	Omnipeak_Peaks
${meta.id}	\$PEAK_COUNT
EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        omnipeak: "1.5"
    END_VERSIONS
    """
}
