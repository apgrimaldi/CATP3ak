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
    path "*.png"                       , emit: png, optional: true
    path "*_mqc.html"                  , emit: mqc_html, optional: true
    path "versions.yml"                , emit: versions

    script:
    // Usiamo shell invece di script per gestire meglio i caratteri speciali di Bash
    shell:
    '''
    #!/usr/bin/env Rscript
    library(DiffBind)

    samples <- read.csv("!{samplesheet}")
    samples$bamReads <- basename(as.character(samples$bamReads))
    samples$Peaks    <- basename(as.character(samples$Peaks))
    if ("bamControl" %in% colnames(samples)) {
        samples$bamControl <- basename(as.character(samples$bamControl))
    }

    db_obj <- dba(sampleSheet=samples)
    
    sample_info <- dba.show(db_obj)
    keep_mask <- as.numeric(sample_info$Intervals) > 0
    
    if(sum(keep_mask) < length(keep_mask)) {
        db_obj <- dba(db_obj, mask=keep_mask)
    }

    pdf("diffbind_correlation.pdf")
    plot(db_obj)
    dev.off()

    png("diffbind_correlation.png", width=800, height=800, res=120)
    plot(db_obj)
    dev.off()

    db_obj <- dba.count(db_obj, bParallel=TRUE)

    pdf("diffbind_profile.pdf")
    try(dba.plotProfile(db_obj, bUseSampleSheet=TRUE))
    dev.off()

    png("diffbind_profile.png", width=1000, height=800, res=120)
    try(dba.plotProfile(db_obj, bUseSampleSheet=TRUE))
    dev.off()

    # Passiamo a Bash per generare gli HTML
    system("
        if [ -f diffbind_correlation.png ]; then
            IMG_CORR=$(base64 -w 0 diffbind_correlation.png)
            echo '<div style=\"text-align:center;\"><img src=\"data:image/png;base64,'$IMG_CORR'\" style=\"max-width:100%;\"></div>' > diffbind_corr_mqc.html
        fi

        if [ -f diffbind_profile.png ]; then
            IMG_PROF=$(base64 -w 0 diffbind_profile.png)
            echo '<div style=\"text-align:center;\"><img src=\"data:image/png;base64,'$IMG_PROF'\" style=\"max-width:100%;\"></div>' > diffbind_profile_mqc.html
        fi
    ")

    writeLines(c(paste0("\\"!{task.process}\\":"), paste0("    diffbind: ", packageVersion("DiffBind"))), "versions.yml")
    '''
}
