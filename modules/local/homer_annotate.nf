process HOMER_ANNOTATEPEAKS {
    tag "$meta.id"
    label 'process_medium'
    container 'quay.io/biocontainers/homer:4.11--pl526hc9558a2_3'

    // Rimuoviamo temporaneamente variabili dal publishDir se ne avevi nel config
    publishDir "${params.outdir}/05_annotation", mode: 'copy'

    input:
    tuple val(meta), path(peak)
    path fasta
    path gtf

    output:
    // Usiamo una wildcard semplice per l'output
    tuple val(meta), path("*.txt"), emit: txt
    path "versions.yml"           , emit: versions

    script:
    // Definiamo il nome del file in modo statico e sicuro
    def sample_id = meta.id
    """
    annotatePeaks.pl \\
        $peak \\
        $fasta \\
        -gtf $gtf \\
        -cpu $task.cpus \\
        > ${sample_id}.annotatePeaks.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        homer: 4.11
    END_VERSIONS
    """
}
