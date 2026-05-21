process HOMER_ANNOTATEPEAKS {
    tag "$meta.id"
    label 'process_medium'
    container 'quay.io/biocontainers/homer:4.11--pl526hc9558a2_3'

    input:
    tuple val(meta), path(peak)
    path  fasta
    path  gtf

    output:
    tuple val(meta), path("*.annotatePeaks.txt") , emit: txt
    tuple val(meta), path("*.homer_stats_mqc.txt"), emit: stats_mqc
    path  "versions.yml"                         , emit: versions

    script:
    def args = task.ext.args ?: ''
    
    def peak_name = (peak != null && peak.name != null) ? peak.name : "unknown"
    def type = peak_name.contains('narrow') ? 'narrow' : 'broad'
    
    def prefix = "${meta.id}.${type}"
    def VERSION = '4.11' 
    """
    # 1. Gestione GTF 
    GTF_FILE=\$(basename ${gtf})
    if [[ \${GTF_FILE} == *.gz ]]; then
        gunzip -c ${gtf} > reference.gtf
        ANNOT_FILE="reference.gtf"
    else
        ANNOT_FILE="${gtf}"
    fi

    # 2. Esecuzione Annotazione
    annotatePeaks.pl \\
        $peak \\
        $fasta \\
        -gtf \${ANNOT_FILE} \\
        $args \\
        -cpu $task.cpus \\
        > ${prefix}.annotatePeaks.txt

    # 3.  Statistiche per MultiQC (Custom Content)
    echo "# id: 'homer_annotations_${type}'" > ${prefix}.homer_stats_mqc.txt
    echo "# section_name: 'HOMER Peak Annotation (${type.capitalize()})'" >> ${prefix}.homer_stats_mqc.txt
    echo "# format: 'tsv'" >> ${prefix}.homer_stats_mqc.txt
    echo "# plot_type: 'bargraph'" >> ${prefix}.homer_stats_mqc.txt
    echo "# pconfig:" >> ${prefix}.homer_stats_mqc.txt
    echo "#    title: 'Peak Annotation Distribution (${type})'" >> ${prefix}.homer_stats_mqc.txt
    echo -e "Sample\\tIntergenic\\tTTS\\texon\\tintron\\tpromoter-TSS" >> ${prefix}.homer_stats_mqc.txt
    
    # Estrazione conteggi
    INTERGENIC=\$(cut -f8 ${prefix}.annotatePeaks.txt | grep -c "Intergenic" || true)
    TTS=\$(cut -f8 ${prefix}.annotatePeaks.txt | grep -c "TTS" || true)
    EXON=\$(cut -f8 ${prefix}.annotatePeaks.txt | grep -c "exon" || true)
    INTRON=\$(cut -f8 ${prefix}.annotatePeaks.txt | grep -c "intron" || true)
    PROMOTER=\$(cut -f8 ${prefix}.annotatePeaks.txt | grep -i -c "promoter-TSS" || true)
    
    echo -e "${prefix}\\t\$INTERGENIC\\t\$TTS\\t\$EXON\\t\$INTRON\\t\$PROMOTER" >> ${prefix}.homer_stats_mqc.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        homer: $VERSION
    END_VERSIONS
    """
}
