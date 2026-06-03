process DIFFBIND {
    tag "diffbind_analysis"
    label 'process_high'
    container 'quay.io/biocontainers/bioconductor-diffbind:3.20.0--r45ha27e39d_0'

    input:
    val prefix      // <--- Parametro aggiunto per differenziare MACS3 da Lanceotron
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

    db_obj <- dba(sampleSheet=samples)
    
    sample_info <- dba.show(db_obj)
    keep_mask <- as.numeric(sample_info\$Intervals) > 0
    if(sum(keep_mask) < length(keep_mask)) {
        db_obj <- dba(db_obj, mask=keep_mask)
    }

    # --- Aggiunto prefisso a tutti i file ---
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

    # 1. CONTEGGIO DEI READ SUI PICCHI (FIX BLACKLIST)
    # Aggiunto tryCatch per evitare che il fallimento di GenomeInfoDb faccia crashare lo script
    tryCatch({
        db_obj <- dba.count(db_obj, bParallel=TRUE, bUseSummarizeOverlaps=FALSE)
    }, error = function(e) {
        print(paste("Errore nel conteggio iniziale, tento il fallback:", e\$message))
        db_obj <<- dba.count(db_obj, bParallel=FALSE, bUseSummarizeOverlaps=FALSE)
    })

    # === ESTRAZIONE E SALVATAGGIO MATRICE DEI CONTEGGI NORMALIZZATI ===
    try({
        counts_data <- dba.peakset(db_obj, bRetrieve=TRUE)
        if (!is.null(counts_data)) {
            write.csv(as.data.frame(counts_data), "${prefix}_diffbind_counts_matrix.csv", row.names=FALSE)
        }
    }, silent=FALSE) 
    # =======================================================================

    try({
        cor_matrix <- dba.overlap(db_obj, mode=DBA_OL_COR)
        write.table(cor_matrix, file="${prefix}_diffbind_correlation_mqc.txt", sep="\t", quote=FALSE, col.names=NA)
    }, silent=TRUE)

    analysis_status <- try({
        contrast_category <- if ("Condition" %in% colnames(samples) && length(unique(samples\$Condition)) > 1) DBA_CONDITION else DBA_ANTIBODY
        
        # [TEST DIAGNOSTICO] Abbassato minMembers a 1 per salvare i campioni con bassa sovrapposizione
        db_obj <- dba.contrast(db_obj, categories=contrast_category, minMembers=1)
        
        # [FIX DESEQ2] Forziamo esplicitamente l'uso di DESeq2
        db_obj <- dba.analyze(db_obj, method=DBA_DESEQ2)
    }, silent=TRUE)

    # Inizializziamo un file vuoto di sicurezza se non ci fossero picchi significativi
    file.create("${prefix}_diffbind_sig_peaks.bed")

    if (!inherits(analysis_status, "try-error") && !is.null(db_obj\$contrasts)) {
        res_db <- dba.report(db_obj)
        
        if(!is.null(res_db)) {
            write.csv(as.data.frame(res_db), "${prefix}_diff_bind_results.csv")

            # ESTRAZIONE E SCRITTURA DELLA MATRICE FILTRATA IN BED PER PROFILEPLYR
            if(length(res_db) > 0) {
                df_sig <- as.data.frame(res_db)
                bed_sig <- df_sig[, c("seqnames", "start", "end")]
                bed_sig\$name <- paste0("DB_site_", 1:nrow(bed_sig))
                bed_sig\$score <- df_sig\$FDR
                write.table(bed_sig, "${prefix}_diffbind_sig_peaks.bed", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
            }
        }

        # Generazione PCA solo se l'analisi differenziale ha prodotto risultati
        try({
            png("${prefix}_diffbind_pca.png", width=1000, height=800, res=150)
            dba.plotPCA(db_obj, attributes=contrast_category, label=DBA_ID)
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

    writeLines(c(
        "\\"${task.process}\\":",
        paste0("    diffbind: ", packageVersion("DiffBind"))
    ), "versions.yml")
    """
}
