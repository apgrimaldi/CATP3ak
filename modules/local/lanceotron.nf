process LANCEOTRON {
    tag "${meta.id}"
    label 'process_high'
    // Utilizziamo il container ufficiale di LanceOtron
    container 'quay.io/biocontainers/lanceotron:1.2.7--pyhdfd78af_0'

    publishDir "${params.outdir}/05_peak_calling/lanceotron", mode: 'copy'

    input:
    // Il modulo ora riceve Meta, BAM IP, BW IP, e opzionalmente BAM CTRL e BW CTRL
    tuple val(meta), path(bam_ip), path(bw_ip), path(bam_ctrl), path(bw_ctrl)

    output:
    tuple val(meta), path("*_peaks.bed")      , emit: peaks
    tuple val(meta), path("*_counts.txt")     , emit: counts_mqc, optional: true
    path "versions.yml"                       , emit: versions

    script:
    def prefix = "${meta.id}"
    
    // Logica per gestire il controllo: se presente (ChIP-seq), lo aggiungiamo al comando
    def control_cmd = (bam_ctrl && bw_ctrl) ? "--control_bam ${bam_ctrl} --control_bigwig ${bw_ctrl}" : ""
    
    """
    # Peak calling con Deep Learning (LanceOtron)
    # Usiamo direttamente i BigWig a 1bp generati da DeepTools
    # 'call_peaks' è il comando standard per l'analisi completa
    
    lanceotron call_peaks ${bam_ip} \\
        --bigwig ${bw_ip} \\
        ${control_cmd} \\
        --output_directory . \\
        --threshold 0.9 \\
        --window_size 1000

    # Rinominiamo gli output per coerenza con il workflow
    # Lanceotron solitamente genera 'L_extract_peaks.bed'
    if [ -f L_extract_peaks.bed ]; then
        mv L_extract_peaks.bed ${prefix}_lanceotron_peaks.bed
    fi

    # Generiamo un file di conteggio semplice per MultiQC
    echo "Sample Peaks" > ${prefix}.lanceotron_counts.txt
    if [ -f ${prefix}_lanceotron_peaks.bed ]; then
        COUNT=\$(grep -v "^#" ${prefix}_lanceotron_peaks.bed | wc -l)
        echo "${prefix} \$COUNT" >> ${prefix}.lanceotron_counts.txt
    else
        echo "${prefix} 0" >> ${prefix}.lanceotron_counts.txt
    fi

    cat <<EOF > versions.yml
    "${task.process}":
        lanceotron: 1.2.7
    EOF
    """
}
