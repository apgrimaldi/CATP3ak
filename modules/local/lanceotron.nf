process LANCEOTRON {
    tag "${meta.id}"
    label 'process_high'
    container 'quay.io/biocontainers/lanceotron:1.2.7--pyhdfd78af_0'

    publishDir "${params.outdir}/05_peak_calling/lanceotron", mode: 'copy'

    input:
    tuple val(meta), path(bam_ip), path(bw_ip), path(bam_ctrl), path(bw_ctrl)

    output:
    tuple val(meta), path("*_peaks.bed")      , emit: peaks
    tuple val(meta), path("*_counts.txt")     , emit: counts_mqc, optional: true
    path "versions.yml"                       , emit: versions

    script:
    def prefix = "${meta.id}"
    
    if (bam_ctrl && bw_ctrl) {
        """
        # ChIP-seq con Controllo
        # file (posizionale) = bw_ip
        # -i = bw_ctrl
        # -f = cartella di output (.)
        lanceotron callPeaksInput \\
            ${bw_ip} \\
            -i ${bw_ctrl} \\
            -f . \\
            -t 0.9 \\
            -w 1000

        # Rinominiamo l'output. Lanceotron genera file basandosi sul nome del BigWig
        # o un generico L_extract_peaks.bed
        if [ -f L_extract_peaks.bed ]; then
            mv L_extract_peaks.bed ${prefix}_lanceotron_peaks.bed
        elif [ -f "${bw_ip.baseName}_peaks.bed" ]; then
            mv "${bw_ip.baseName}_peaks.bed" ${prefix}_lanceotron_peaks.bed
        fi

        # Generazione conteggio per MultiQC
        echo "Sample Peaks" > ${prefix}.lanceotron_counts.txt
        if [ -f ${prefix}_lanceotron_peaks.bed ]; then
            COUNT=\$(grep -v "^#" ${prefix}_lanceotron_peaks.bed | wc -l)
            echo "${prefix} \$COUNT" >> ${prefix}.lanceotron_counts.txt
        else
            echo "${prefix} 0" >> ${prefix}.lanceotron_counts.txt
        fi

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            lanceotron: 1.2.7
        END_VERSIONS
        """
    } else {
        """
        # ATAC-seq (Senza controllo)
        # file (posizionale) = bw_ip
        lanceotron callPeaks \\
            ${bw_ip} \\
            -f . \\
            -t 0.9 \\
            -w 1000

        if [ -f L_extract_peaks.bed ]; then
            mv L_extract_peaks.bed ${prefix}_lanceotron_peaks.bed
        elif [ -f "${bw_ip.baseName}_peaks.bed" ]; then
            mv "${bw_ip.baseName}_peaks.bed" ${prefix}_lanceotron_peaks.bed
        fi

        echo "Sample Peaks" > ${prefix}.lanceotron_counts.txt
        if [ -f ${prefix}_lanceotron_peaks.bed ]; then
            COUNT=\$(grep -v "^#" ${prefix}_lanceotron_peaks.bed | wc -l)
            echo "${prefix} \$COUNT" >> ${prefix}.lanceotron_counts.txt
        else
            echo "${prefix} 0" >> ${prefix}.lanceotron_counts.txt
        fi

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            lanceotron: 1.2.7
        END_VERSIONS
        """
    }
}
