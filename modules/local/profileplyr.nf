process PROFILEPLYR {
    tag "${label}"
    label 'process_high'
    container 'quay.io/biocontainers/bioconductor-profileplyr:1.22.0--r44hdfd78af_0'

    publishDir "results/profile_heatmaps", mode: 'copy'

    input:
    path(diff_peaks, stageAs: 'diff_peaks/*') // [INPUT 1] I picchi di DiffBind
    path(raw_peaks,  stageAs: 'raw_peaks/*')  // [INPUT 2] I picchi grezzi (il piano B)
    path(bigwigs,    stageAs: 'bigwigs/*')    // [INPUT 3] I BigWig
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

    # IL FIX DEFINITIVO PER IL CRASH "CAIRO": 
    # Impedisce a R di tentare la compressione dell'immagine quando i picchi sono migliaia.
    ht_opt(use_raster = FALSE)

    print(paste("Analisi Profileplyr:", "${label}"))

    diff_files <- list.files("diff_peaks", pattern = "\\\\.(bed|narrowPeak|broadPeak)\$", full.names = TRUE)
    raw_files  <- list.files("raw_peaks", pattern = "\\\\.(bed|narrowPeak|broadPeak)\$", full.names = TRUE)
    bw_files   <- list.files("bigwigs", pattern = "\\\\.(bw|bigWig)\$", full.names = TRUE)

    # Funzione per caricare e unire i picchi
    import_and_merge_peaks <- function(files) {
        if (length(files) == 0) return(NULL)
        gr_list <- lapply(files, function(f) {
            tryCatch({
                gr <- rtracklayer::import(f)
                if (!is.null(mcols(gr)\$score)) mcols(gr)\$score <- as.numeric(mcols(gr)\$score)
                return(gr)
            }, error = function(e) NULL)
        })
        gr_list <- gr_list[!sapply(gr_list, is.null)]
        if (length(gr_list) == 0) return(NULL)
        all_peaks <- do.call(c, gr_list)
        return(all_peaks)
    }

    # ---- LOGICA FALLBACK AUTOMATICO (COME PRIMA) ----
    peaks_gr <- import_and_merge_peaks(diff_files)
    used_mode <- "Picchi Differenziali (DiffBind)"

    if (is.null(peaks_gr) || length(peaks_gr) < 3) {
        print("Picchi differenziali insufficienti o assenti. Attivo il FALLBACK sui picchi totali/grezzi.")
        peaks_gr <- import_and_merge_peaks(raw_files)
        used_mode <- "Picchi Totali/Grezzi (Fallback)"
    } else {
        print("Picchi differenziali trovati! Utilizzo i risultati di DiffBind.")
    }
    # -------------------------------------------------

    if (is.null(peaks_gr) || length(peaks_gr) < 3) {
        cat(paste0(
            "\\n",
            "<div style='text-align: center; padding: 20px; border: 1px solid #cc0000;'>\\n",
            "  <h3>Profile Analysis: ", "${label}", "</h3>\\n",
            "  <p style='color: #cc0000;'>Nessun picco valido trovato come input. Heatmap non generata.</p>\\n",
            "</div>"
        ), file="${label}_profileplyr_mqc.html")
        
        file.create("${label}_profile_heatmap.png")
        file.create("${label}_profile_heatmap.pdf")
    } else {
        
        # ORDINAMENTO PER SCORE
        if (!is.null(mcols(peaks_gr)\$score)) {
            peaks_gr <- peaks_gr[order(mcols(peaks_gr)\$score, decreasing = TRUE)]
        }

        # FIX RAM: Limita a 15.000 picchi per non far esplodere la memoria
        if (length(peaks_gr) > 15000) {
            peaks_gr <- head(peaks_gr, 15000)
        }

        # Fusione picchi sovrapposti
        peaks_gr <- GenomicRanges::reduce(peaks_gr)

        # 2. Creazione oggetto Profileplyr
        rtracklayer::export(peaks_gr, "merged_testRanges.bed")

        pro_chip <- BamBigwig_to_chipProfile(
            signalFiles = bw_files,
            testRanges = "merged_testRanges.bed",
            format = "bigwig",
            style = "point",
            bin_size = 50,
            distanceAround = 2000
        )

        pro_obj <- as_profileplyr(pro_chip)
        sampleData(pro_obj)\$sample_id <- sub("\\\\.(bw|bigWig)\$", "", basename(bw_files))

        # 3. Plotting sicuro
        tryCatch({
            ht <- generateEnrichedHeatmap(pro_obj)
            
            # PDF ad alta risoluzione
            pdf("${label}_profile_heatmap.pdf", width=8, height=10)
            print(ht)
            dev.off()

            # PNG per MultiQC (con doppio tentativo per massima stabilità)
            png_success <- FALSE
            tryCatch({
                png("${label}_profile_heatmap.png", width=1200, height=1400, res=150, type="cairo")
                print(ht)
                dev.off()
                png_success <- TRUE
            }, error = function(e_cairo) {
                tryCatch({
                    png("${label}_profile_heatmap.png", width=1200, height=1400, res=150)
                    print(ht)
                    dev.off()
                    png_success <- TRUE
                }, error = function(e_nocairo) {
                    png_success <- FALSE
                })
            })

            if (png_success) {
                img_64 <- base64encode("${label}_profile_heatmap.png")
                cat(paste0(
                    "\\n",
                    "<div style='text-align: center; padding: 20px;'>\\n",
                    "  <h3>Profile Analysis (", "${label}", ")</h3>\\n",
                    "  <p><strong>Modalità usata: </strong> ", used_mode, "</p>\\n",
                    "  <img src='data:image/png;base64,", img_64, "' style='width: 600px; max-width: 100%; height: auto; border: 1px solid #ddd;'>\\n",
                    "</div>"
                ), file="${label}_profileplyr_mqc.html")
            } else {
                stop("Impossibile generare il PNG.")
            }

        }, error = function(e) {
            cat(paste0(
