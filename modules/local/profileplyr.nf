process PROFILEPLYR {
    tag "${label}"
    label 'process_high'
    container 'quay.io/biocontainers/bioconductor-profileplyr:1.22.0--r44hdfd78af_0'

    publishDir "results/profile_heatmaps", mode: 'copy'

    input:
    path(diff_peaks, stageAs: 'diff_peaks/*')
    path(raw_peaks,  stageAs: 'raw_peaks/*')
    path(bigwigs,    stageAs: 'bigwigs/*')
    val label

    output:
    path "*_profile_heatmap.pdf"       , emit: pdf, optional: true
    path "*_profile_heatmap.png"       , emit: png, optional: true
    path "*_profileplyr_mqc.html"      , emit: mqc_html, optional: true
    path "versions.yml"                , emit: versions

    script:
    """
    #!/usr/bin/env Rscript
    library(profileplyr)
    library(base64enc)
    library(rtracklayer)
    library(GenomicRanges)
    library(ComplexHeatmap)

    print(paste("Analisi Profileplyr:", "${label}"))

    # Legge i file dalle cartelle di staging
    diff_files <- list.files("diff_peaks", pattern = "\\\\.(bed|narrowPeak|broadPeak)\$", full.names = TRUE)
    raw_files  <- list.files("raw_peaks", pattern = "\\\\.(bed|narrowPeak|broadPeak)\$", full.names = TRUE)
    bw_files   <- list.files("bigwigs", pattern = "\\\\.(bw|bigWig)\$", full.names = TRUE)

    # Funzione sicura per importare e unire i picchi
    import_and_merge_peaks <- function(files) {
        if (length(files) == 0) return(NULL)
        gr_list <- lapply(files, function(f) {
            tryCatch({
                # Controlla che il file non sia vuoto (size > 0)
                if (file.info(f)\$size == 0) return(NULL)
                
                gr <- rtracklayer::import(f)
                if (!is.null(mcols(gr)\$score)) mcols(gr)\$score <- as.numeric(mcols(gr)\$score)
                return(gr)
            }, error = function(e) NULL)
        })
        gr_list <- gr_list[!sapply(gr_list, is.null)]
        if (length(gr_list) == 0) return(NULL)
        return(do.call(c, gr_list))
    }

    # 1. TENTA CON I PICCHI DIFFBIND
    peaks_gr <- import_and_merge_peaks(diff_files)
    used_mode <- "Picchi Differenziali (DiffBind)"

    # 2. FALLBACK SUI PICCHI ESTERNI (MACS3/Lanceotron)
    if (is.null(peaks_gr) || length(peaks_gr) < 3) {
        peaks_gr <- import_and_merge_peaks(raw_files)
        used_mode <- "Picchi Totali/Grezzi (Fallback)"
    } 

    # 3. SE FALLISCE ANCHE IL FALLBACK -> ESCE IN MODO PULITO
    if (is.null(peaks_gr) || length(peaks_gr) < 3) {
        cat(paste0(
            "\\n<div style='text-align: center; padding: 20px; border: 1px solid #cc0000;'>\\n",
            "  <h3>Profile Analysis: ", "${label}", "</h3>\\n",
            "  <p style='color: #cc0000;'>Nessun picco valido trovato come input. Heatmap saltata.</p>\\n",
            "</div>"
        ), file="${label}_profileplyr_mqc.html")
        file.create("${label}_profile_heatmap.png")
        file.create("${label}_profile_heatmap.pdf")
        quit(save="no", status=0) # Esce senza errore per Nextflow
    }

    # === BLOCCO TRY-CATCH GLOBALE (Niente può far crashare la pipeline da qui in poi) ===
    tryCatch({
        
        # Ordina per score (prendendo il valore assoluto per non perdere i Loss)
        if (!is.null(mcols(peaks_gr)\$score)) {
            peaks_gr <- peaks_gr[order(abs(mcols(peaks_gr)\$score), decreasing = TRUE)]
        }

        # Sicurezza memoria: max 10.000 picchi
        if (length(peaks_gr) > 10000) {
            peaks_gr <- head(peaks_gr, 10000)
        }

        # Logica Gain/Loss solo se veniamo da DiffBind
        if (used_mode == "Picchi Differenziali (DiffBind)" && !is.null(mcols(peaks_gr)\$score)) {
            peak_status <- ifelse(mcols(peaks_gr)\$score > 0, "Gain", "Loss")
            test_ranges <- split(peaks_gr, peak_status)
        } else {
            peaks_gr <- GenomicRanges::reduce(peaks_gr)
            test_ranges <- peaks_gr
        }

        # Calcolo pesante (ora protetto)
        pro_chip <- BamBigwig_to_chipProfile(
            signalFiles = bw_files,
            testRanges = test_ranges,
            format = "bigwig",
            style = "point",
            bin_size = 50,
            distanceAround = 2000
        )

        pro_obj <- as_profileplyr(pro_chip)
        sampleData(pro_obj)\$sample_id <- sub("\\\\.(bw|bigWig)\$", "", basename(bw_files))

        # Generazione Plot
        ht <- generateEnrichedHeatmap(pro_obj, use_raster = FALSE, column_names_gp = grid::gpar(fontsize = 8))
        
        pdf("${label}_profile_heatmap.pdf", width=10, height=12)
        print(ht)
        dev.off()

        png_success <- FALSE
        tryCatch({
            png("${label}_profile_heatmap.png", width=1500, height=1800, res=150, type="cairo")
            print(ht)
            dev.off()
            png_success <- TRUE
        }, error = function(e_cairo) {
            tryCatch({
                png("${label}_profile_heatmap.png", width=1500, height=1800, res=150)
                print(ht)
                dev.off()
                png_success <- TRUE
            }, error = function(e_nocairo) { png_success <- FALSE })
        })

        if (png_success) {
            img_64 <- base64encode("${label}_profile_heatmap.png")
            cat(paste0(
                "\\n<div style='text-align: center; padding: 20px;'>\\n",
                "  <h3>Profile Analysis (", "${label}", ")</h3>\\n",
                "  <p><strong>Modalità usata: </strong> ", used_mode, "</p>\\n",
                "  <img src='data:image/png;base64,", img_64, "' style='width: 800px; max-width: 100%; height: auto; border: 1px solid #ddd;'>\\n",
                "</div>"
            ), file="${label}_profileplyr_mqc.html")
        } else {
            stop("Rendering PNG fallito.")
        }

    }, error = function(e) {
        # SE QUALSIASI COSA FALLISCE (Memoria, Calcolo, Plot), CATTURA L'ERRORE E SALVA LA PIPELINE
        cat(paste0(
            "\\n<div style='text-align: center; padding: 20px; border: 1px solid #ffcc00;'>\\n",
            "  <h3>Profile Analysis: ", "${label}", "</h3>\\n",
            "  <p><strong>Avviso:</strong> Il grafico non è stato generato a causa di un problema tecnico (es. memoria insufficiente o dati non idonei).</p>\\n",
            "  <p><em>Dettaglio tecnico: ", e\$message, "</em></p>\\n",
            "</div>"
        ), file="${label}_profileplyr_mqc.html")
        file.create("${label}_profile_heatmap.png")
        file.create("${label}_profile_heatmap.pdf")
    })

    writeLines(c("\\"${task.process}\\":", paste0("    profileplyr: ", packageVersion("profileplyr"))), "versions.yml")
    """
}
