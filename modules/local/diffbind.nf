process DIFFBIND {
    tag "diffbind_analysis_lanceotron"
    label 'process_high'
    container 'quay.io/biocontainers/bioconductor-diffbind:3.20.0--r45ha27e39d_0'

    input:
    val prefix
    path samplesheet
    path bams
    path bais
    path peaks

    output:
    path "*.pdf"                       , emit: pdf, optional: true
    path "*.csv"                       , emit: csv, optional: true
    path "*_mqc.html"                  , emit: mqc_html, optional: true
    path "*_correlation_mqc.txt"       , emit: mqc_txt, optional: true
    path "*.png"                       , emit: png, optional: true
    path "*_sig_peaks.bed"             , emit: sig_bed, optional: true 
    path "versions.yml"                , emit: versions

    script:
    """
    #!/usr/bin/env Rscript
    library(DiffBind)
    library(base64enc)
    library(BiocParallel)

    allocated_cores <- as.numeric("${task.cpus}")
    register(MulticoreParam(workers = allocated_cores, progressbar = TRUE))
    options(cores = allocated_cores)

    samples <- read.csv("${samplesheet}")
    samples\$bamReads <- basename(as.character(samples\$bamReads))
    samples\$Peaks    <- basename(as.character(samples\$Peaks))
    if ("bamControl" %in% colnames(samples)) {
        samples\$bamControl <- basename(as.character(samples\$bamControl))
    }

    if ("antibody" %in% colnames(samples)) {
        colnames(samples)[colnames(samples) == "antibody"] <- "Factor"
    }

    db_obj <- dba(sampleSheet=samples)
    
    sample_info <- dba.show(db_obj)
    keep_mask <- as.numeric(sample_info\$Intervals) > 0
    if(sum(keep_mask) < length(keep_mask)) {
        db_obj <- dba(db_obj, mask=keep_mask)
    }

    png("${prefix}_diffbind_correlation.png", width=1000, height=1000, res=150)
    plot(db_obj, margin=20)
    dev.off()

    img_corr_64 <- base64encode("${prefix}_diffbind_correlation.png")
    
    cat(paste0(
        "\\n",
        "<div style='text-align: center; padding: 20px;'>\\n",
        "  <img src='data:image/png;base64,", img_corr_64, "' style='width: 500px; max-width: 100%; height: auto;'>\\n",
        "</div>"
    ), file="${prefix}_diffbind_corr_mqc.html")

    db_obj <- dba.blacklist(db_obj, blacklist=FALSE, greylist=FALSE)
    db_obj <- dba.count(db_obj, bParallel=TRUE, bUseSummarizeOverlaps=FALSE)

    try({
        counts_data <- dba.peakset(db_obj, bRetrieve=TRUE)
        if (!is.null(counts_data)) {
            write.csv(as.data.frame(counts_data), "${prefix}_diffbind_counts_matrix.csv", row.names=FALSE)
        }
    }, silent=TRUE) 

    try({
        cor_matrix <- dba.overlap(db_obj, mode=DBA_OL_COR)
        write.table(cor_matrix, file="${prefix}_diffbind_correlation_mqc.txt", sep="\t", quote=FALSE, col.names=NA)
    }, silent=TRUE)

    analysis_status <- tryCatch({
        db_tmp <- dba.contrast(db_obj, categories=DBA_FACTOR, minMembers=2)
        dba.analyze(db_tmp, method=DBA_DESEQ2)
    }, error = function(e1) {
        tryCatch({
            db_tmp <- dba.contrast(db_obj, categories=DBA_FACTOR, minMembers=1)
            dba.analyze(db_tmp, method=DBA_DESEQ2)
        }, error = function(e2) {
            return(NULL)
        })
    })

    file.create("${prefix}_diffbind_sig_peaks.bed")

    if (!is.null(analysis_status) && !is.null(analysis_status\$contrasts)) {
        res_db <- dba.report(analysis_status)
        
        if(!is.null(res_db)) {
            write.csv(as.data.frame(res_db), "${prefix}_diff_bind_results.csv")

            if(length(res_db) > 0) {
                df_sig <- as.data.frame(res_db)
                bed_sig <- df_sig[, c("seqnames", "start", "end")]
                bed_sig\$name <- paste0("DB_site_", 1:nrow(bed_sig))
                bed_sig\$score <- df_sig\$FDR
                write.table(bed_sig, "${prefix}_diffbind_sig_peaks.bed", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
            }

            try({
                png("${prefix}_diffbind_pca.png", width=1000, height=800, res=150)
                dba.plotPCA(analysis_status, attributes=DBA_FACTOR, label=DBA_ID)
                dev.off()

                img_pca_64 <- base64encode("${prefix}_diffbind_pca.png")
                
                cat(paste0(
                    "\\n",
                    "<div style='text-align: center; padding: 20px;'>\\n",
                    "  <img src='data:image/png;base64,", img_pca_64, "' style='width: 500px; max-width: 100%; height: auto;'>\\n",
                    "</div>"
                ), file="${prefix}_diffbind_pca_mqc.html")
            }, silent=TRUE)
        }
    }

    writeLines(c(
        "\\"${task.process}\\":",
        paste0("    diffbind: ", packageVersion("DiffBind"))
    ), "versions.yml")
    """
}
