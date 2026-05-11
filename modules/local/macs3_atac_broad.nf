process MACS3_ATAC_BROAD {
    tag "$meta.id"
    label 'process_medium'
    container 'quay.io/biocontainers/macs3:3.0.1--py311h0152c62_3'

    input:
    tuple val(meta), path(bam)
    val gsize // Riceve il valore corretto dal workflow (automatico da iGenomes o manuale)

    output:
    tuple val(meta), path("*.broadPeak") , emit: peaks
    tuple val(meta), path("*.broad_counts.txt"), emit: count_broad // Cambiato in tuple per coerenza con gli altri moduli
    path "versions.yml"                  , emit: versions

    script:
    def prefix   = "${meta.id}_atac_broad"
    def format   = meta.single_end ? 'BAM' : 'BAMPE'

    """
    macs3 callpeak \\
        -t $bam \\
        -f $format \\
        -g $gsize \\
        -n $prefix \\
        --nomodel --shift -100 --extsize 200 \\
        --broad \\
        --broad-cutoff 0.1

    # Estrazione automatica del conteggio per il grafico MultiQC
    if [ -f ${prefix}_peaks.broadPeak ]; then
        count=\$(wc -l < ${prefix}_peaks.broadPeak)
    else
        count=0
    fi
    
    # Generazione file per MultiQC con header per chiarezza
    echo -e "Sample\\tBroad_Peaks" > ${prefix}.broad_counts.txt
    echo -e "${meta.id}\\t\$count" >> ${prefix}.broad_counts.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        macs3: \$(macs3 --version | sed 's/macs3 //g')
    END_VERSIONS
    """
}
