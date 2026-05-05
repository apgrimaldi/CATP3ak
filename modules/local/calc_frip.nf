process CALC_FRIP {
    tag "$meta.id"
    label 'process_medium'

    container 'quay.io/biocontainers/subread:2.0.1--hed695b0_0'

    input:
    tuple val(meta), path(bam), path(peak)

    output:
    tuple val(meta), path("*.summary"), emit: summary
    path "*.frip.txt"                 , emit: txt
    path "versions.yml"               , emit: versions

    script:
    def prefix = "${meta.id}"
    """
    # Trasformiamo il file peak BED in un formato SAF (richiesto da featureCounts)
    # SAF: GeneID  Chr  Start  End  Strand
    awk 'BEGIN{FS=OFS="\t"; print "GeneID\tChr\tStart\tEnd\tStrand"} {print "peak_"NR, \$1, \$2+1, \$3, "."}' $peak > peaks.saf

    featureCounts \\
        -p \\
        -a peaks.saf \\
        -F SAF \\
        -o ${prefix}.featureCounts.txt \\
        $bam

    # Calcoliamo la FRiP manualmente per un file di testo veloce
    READS_IN_PEAKS=\$(grep -v "Status" ${prefix}.featureCounts.txt.summary | awk '{sum+=\$2} END {print sum}')
    TOTAL_READS=\$(samtools view -c $bam) # Nota: richiede samtools nel container o un pipe
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        subread: 2.0.1
    END_VERSIONS
    """
}
