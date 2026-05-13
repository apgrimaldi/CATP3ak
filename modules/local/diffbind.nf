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
    path "*.pdf"                       , emit: pdf, optional: true
    path "*.csv"                       , emit: csv, optional: true
    path "*_mqc.html"                  , emit: mqc_html, optional: true
    path "diffbind_correlation_mqc.txt", emit: mqc_txt, optional: true
    path "*.png"                       , emit: png, optional: true
    path "versions.yml"                , emit: versions

    script:
    """
    #!/usr/bin/env Rscript
    library(DiffBind)

    samples <- read.csv("${samplesheet}")
    samples\$bamReads <- basename(as.character(samples\$bamReads))
    samples\$Peaks    <- basename(as.character(samples\$Peaks))
    if ("bamControl" %in% colnames(samples)) {
        samples\$bamControl <- basename(as.character(samples\$bamControl))
    }

    db_obj <- dba(sampleSheet=samples)
    
    pdf("diffbind_correlation.pdf")
    plot(db_obj)
    dev.off()

    png("diffbind_correlation.png", width=800, height=800, res=120)
    plot(db_obj)
    dev.off()

    cat(paste0(
        "\\n",
        "<div style='text-align: center;'>\\n",
        "  <img src='diffbind_correlation.png' style='max-width: 100%; height: auto;'>\\n",
        "</div>"
    ), file="diffbind_corr_mqc.html")

    db_obj <- dba.count(db_obj, bParallel=TRUE)

    try({
        cor_matrix <- dba.overlap(db_obj, mode=DBA_OL_COR)
        write.table(cor_matrix, file="diffbind_correlation_mqc.txt", sep="\t", quote=FALSE, col.names=NA)
    }, silent=TRUE)

    analysis_status <- try({
        contrast_category <- if ("Condition" %in% colnames(samples) && length(unique(samples\$Condition)) > 1) DBA_CONDITION else DBA_ANTIBODY
        db_obj <- dba.contrast(db_obj, categories=contrast_category, minMembers=2)
        db_obj <- dba.analyze(db_obj)
    }, silent=FALSE)

    if (!inherits(analysis_status, "try-error") && !is.null(db_obj\$contrasts)) {
        res_db <- dba.report(db_obj)
        write.csv(as.data.frame(res_db), "diff_bind_results.csv")

        png("diffbind_pca.png", width=1000, height=800, res=120)
        dba.plotPCA(db_obj, attributes=contrast_category, label=DBA_ID)
        dev.off()

        pdf("diffbind_volcano.pdf")
        dba.plotVolcano(db_obj)
        dev.off()
    }

    writeLines(c(
        "\\"${task.process}\\":",
        paste0("    diffbind: ", packageVersion("DiffBind"))
    ), "versions.yml")
    """
}
