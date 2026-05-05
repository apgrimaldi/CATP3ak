process MULTIQC {
    label 'process_medium'
    container 'quay.io/biocontainers/multiqc:1.21--pyhdfd78af_0'

    publishDir "${params.outdir}/00_MultiQC", mode: 'copy'

    input:
    path multiqc_config      // 1
    path workflow_summary    // 2
    path ('fastqc/*')        // 3
    path ('trimgalore/*')    // 4
    path ('alignment/*')     // 5
    path ('picard/*')        // 6
    path ('samtools/*')      // 7
    path ('deeptools/*')     // 8  <-- Assicurati che questo ci sia
    path ('macs3/*')         // 9
    path ('frip/*')          // 10
    path ('homer/*')         // 11 <-- Aggiunto/Verificato
    path versions            // 12 <-- Questo è il dodicesimo

    output:
    path "*multiqc_report.html", emit: report
    path "*_data"              , emit: data
    path "versions.yml"        , emit: versions

    script:
    def args = task.ext.args ?: ''
    def config = multiqc_config.name != 'empty_config' ? "--config $multiqc_config" : ''
    """
    multiqc \\
        -f \\
        $args \\
        $config \\
        .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$(multiqc --version | sed 's/multiqc, version //g')
    END_VERSIONS
    """
}
