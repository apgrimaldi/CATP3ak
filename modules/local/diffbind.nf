process DIFFBIND {
    tag "diffbind_analysis"
    label 'process_high'
    container 'quay.io/biocontainers/bioconductor-diffbind:3.20.0--r45ha27e39d_0'

    input:
    path samplesheet
    path bams
    path bais
    path peaks

    output:
    path "*.pdf"                   , emit: pdf, optional: true
    path "*.csv"                   , emit: csv, optional: true
    path "*_mqc.html"              , emit: mqc_html, optional: true
    path "*.png"                   , emit: png, optional: true
    path "versions.yml"            , emit: versions

    script:
    """
    #!/usr/bin/env Rscript
    library(DiffBind)

    # Caricamento campioni
    samples <- read.csv("${samplesheet}")
    samples\$bamReads <- basename(samples\$bamReads)
    samples\$Peaks    <- basename(samples\$Peaks)
    
    if ("ControlID" %in% colnames(samples)) {
        samples\$bamControl <- basename(samples\$bamControl)
    }

    # Creazione oggetto DBA
    db_obj <- dba(sampleSheet=samples)
    
    # Conteggio letture
    db_obj <- dba.count(db_obj, bParallel=TRUE)

    # Tenta di stabilire il contrasto e analizzare
    # Usiamo un blocco try() per evitare che il crash blocchi l'intera pipeline 
    # se i campioni sono insufficienti per l'analisi statistica
    analysis_status <- try({
        db_obj <- dba.contrast(db_obj, categories=DBA_CONDITION, minMembers=2)
        db_obj <- dba.analyze(db_obj)
    }, silent=TRUE)

    # Se l'analisi è riuscita, genera i report
    if (!inherits(analysis_status, "try-error")) {
        res_db <- dba.report(db_obj)
        write.csv(as.data.frame(res_db), "diff_bind_results.csv")

        png("diffbind_profile.png", width=1000, height=800, res=120)
        try(dba.plotProfile(db_obj, contrast=1))
        dev.off()

        pdf("diffbind_profile.pdf", width=10, height=8)
        try(dba.plotProfile(db_obj, contrast=1))
        dev.off()

        cat(paste0(
            "# id: 'diffbind_profile'\\n",
            "# section_name: 'Differential Binding Profile'\\n",
            "<div style='text-align: center;'>\\n",
            "  <img src='diffbind_profile.png' style='max-width: 100%; height: auto;'>\\n",
            "</div>"
        ), file="diffbind_profile_mqc.html")
    } else {
        # Se fallisce, crea un file log invece di crashare
        writeLines("Analisi differenziale non riuscita: campioni insufficienti o nessun contrasto trovato.", "diffbind_warning.log")
    }

    # Versioning
    writeLines(c(
        "\\"${task.process}\\":",
        paste0("    diffbind: ", packageVersion("DiffBind"))
    ), "versions.yml")
    """
}
