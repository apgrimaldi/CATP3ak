process OMNIPEAK {
    tag "$meta.id"
    label 'process_high'
    container 'biohaz/omnipeak:latest'

    input:
    // Sostituito i due path separati con una lista dinamica di BAM
    tuple val(meta), path(bams)
    path chrom_sizes

    output:
    tuple val(meta), path("*_peaks.bed"), emit: peaks
    path "*_omnipeak_mqc.txt"           , emit: counts_mqc
    path "versions.yml"                 , emit: versions

    script:
    // LOGICA INTELLIGENTE: Il primo file è sempre il treatment. 
    // Se c'è un secondo file, allora è il control.
    def ip_bam = bams[0]
    def ctrl_arg = bams.size() > 1 ? "--control \"${bams[1]}\"" : ""

    """
    java --add-modules=jdk.incubator.vector -Xmx16G -jar /home/omnipeak/build/libs/omnipeak-1.5.build.jar analyze \\
        --threads 4 \\
        --treatment "${ip_bam}" \\
        ${ctrl_arg} \\
        --cs ${chrom_sizes} \\
        --peaks "${meta.id}_peaks.bed" \\
        > omnipeak.log 2>&1

    
    if [ -f "${meta.id}_peaks.bed" ]; then
        PEAK_COUNT=\$(wc -l < ${meta.id}_peaks.bed)
    else
        PEAK_COUNT=0
    fi

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
