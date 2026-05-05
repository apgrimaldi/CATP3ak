process CALC_FRIP {
    tag "$meta.id"
    label 'process_medium'

    // Usiamo lo stesso container di nf-core che contiene bedtools e samtools
    container "quay.io/biocontainers/mulled-v2-8186960447c5cb2faa697666dc1e6d919ad23f3e:3127fcae6b6bdaf8181e21a26ae61231030a9fcb-0"

    input:
    tuple val(meta), path(bam), path(peak)

    output:
    tuple val(meta), path("*.FRiP.txt"), emit: frip
    path "versions.yml"                , emit: versions

    script:
    def prefix = "${meta.id}"
    """
    # 1. Conta le reads che intersecano i picchi
    # -bed converte l'output in formato BED (molto leggero)
    # -c conta quante reads cadono in ogni picco
    READS_IN_PEAKS=\$(bedtools intersect -a $peak -b $bam -c | awk -F '\\t' '{sum += \$NF} END {print sum}')

    # 2. Ottieni il numero totale di reads mappate
    samtools flagstat $bam > ${prefix}.flagstat
    TOTAL_MAPPED=\$(grep 'mapped (' ${prefix}.flagstat | grep -v "primary" | head -n 1 | awk '{print \$1}')

    # 3. Calcola lo score FRiP
    awk -v rip="\$READS_IN_PEAKS" -v total="\$TOTAL_MAPPED" -v samp="${prefix}" \\
        'BEGIN {OFS="\\t"; print "Sample", "ReadsInPeaks", "TotalMapped", "FRiP"; \\
        if (total > 0) print samp, rip, total, rip/total; else print samp, rip, total, 0}' > ${prefix}.FRiP.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
        samtools: \$(samtools --version | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """
}
