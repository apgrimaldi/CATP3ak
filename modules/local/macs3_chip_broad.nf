process MACS3_CHIP_BROAD {
    tag "$meta.id"
    label 'process_medium'
    container 'quay.io/biocontainers/macs3:3.0.1--py311h0152c62_3'

    input:
    tuple val(meta), path(ip_bam), path(control_bam)
    val gsize 

    output:
    // Usiamo il glob * per catturare il file corretto indipendentemente dal nome esatto
    tuple val(meta), path("*.broadPeak")       , emit: peaks
    tuple val(meta), path("*.xls")              , emit: xls
    tuple val(meta), path("*.broad_counts.txt"), emit: count_broad
    path "versions.yml"                        , emit: versions

    script:
    def prefix   = "${meta.id}_broad"
    def format   = meta.single_end ? 'BAM' : 'BAMPE'
    
    // Gestione robusta del controllo: se control_bam esiste nel tuple, lo aggiungiamo
    def args_control = control_bam ? "-c ${control_bam}" : ""

    """
    macs3 callpeak \\
        -t $ip_bam \\
        $args_control \\
        -f $format \\
        -g $gsize \\
        -n $prefix \\
        --broad \\
        --broad-cutoff 0.1

    # Conteggio sicuro (MACS3 aggiunge spesso _peaks al nome fornito con -n)
    # Cerchiamo il file .broadPeak indipendentemente dal nome esatto con protezione 2>/dev/null
    PEAK_FILE=\$(ls *.broadPeak 2>/dev/null | head -n 1)
    if [ -n "\$PEAK_FILE" ] && [ -f "\$PEAK_FILE" ]; then
        count=\$(wc -l < "\$PEAK_FILE")
    else
        count=0
    fi

    # Output per MultiQC
    echo -e "Sample\\tBroad_Peaks" > ${prefix}.broad_counts.txt
    echo -e "${meta.id}\\t\$count" >> ${prefix}.broad_counts.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        macs3: \$(macs3 --version | sed 's/macs3 //g')
    END_VERSIONS
    """
}
