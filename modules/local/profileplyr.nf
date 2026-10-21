process PROFILEPLYR {
    tag "${label}"
    label 'process_high'
    container 'quay.io/biocontainers/bioconductor-profileplyr:1.22.0--r44hdfd78af_0'

    publishDir "${params.outdir}/12_profileplyr", mode: 'copy'

    input:
    path(diff_peaks, stageAs: 'diff_peaks/*')
    path(raw_peaks,  stageAs: 'raw_peaks/*')
    path(bigwigs,    stageAs: 'bigwigs/*')
    path(gtf)
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

    print(paste("Analisi Profileplyr (TSS mode):", "${label}"))

    bw_files <- list.files("bigwigs", pattern = "\\\\.(bw|bigWig)\$", full.names = TRUE)

    # 1. Carica il GTF ed estrai i TSS
    gtf_data <- tryCatch({
        import("${gtf}")
    }, error = function(e) { NULL })

    used_mode <- "Transcription Start Sites (TSS)"

    if (is.null(gtf_data) || length(bw_files) == 0) {
        cat(paste0(
            "\\n<div style='text-align: center; padding: 20px; border: 1px solid #cc0000;'>\\n",
            "  <h3>Profile Analysis: ", "${label}", "</h3>\\n",
            "  <p style='color: #cc0000;'>Errore: GTF o BigWig mancanti. Heatmap non generata.</p>\\n",
            "</div>"
        ), file="${label}_profileplyr_mqc.html")
        file.create("${label}_profile_heatmap.png")
        file.create("${label}_profile_heatmap.pdf")
    } else {
        
        # Filtra i trascritti e ottieni il TSS esatto
        transcripts <- gtf_data[gtf_data\$type == "transcript" | gtf_data\$type == "mRNA"]
        if (length(transcripts) == 0) {
            transcripts <- gtf_data[gtf_data\$type == "gene"]
        }
        
        tss_gr <- resize(transcripts, width = 1, fix = "start")
        
        # --- LA FIX È QUI ---
        # Puliamo i metadati del GTF che causano il crash del DataFrame
        mcols(tss_gr) <- NULL 
        tss_gr\$name <- paste0("TSS_", seq_along(tss_gr))
        # --------------------
        
        tss_gr <- unique(tss_gr)

        # Limita ai 15.000 TSS casuali per evitare crash di memoria
        if (length(tss_gr) > 15000) {
            tss_gr <- tss_gr[sample(seq_along(tss_gr), 15000)]
        }

        export(tss_gr, "tss_ranges.bed")

        # 2. Genera il profilo usando i TSS puliti
        tryCatch({
            pro_chip <- BamBigwig_to_chipProfile(
                signalFiles = bw_files,
                testRanges = "tss_ranges.bed",
                format = "bigwig",
                style = "point",
                bin_size = 50,
                distanceAround = 2000
            )

            pro_obj <- as_profileplyr(pro_chip)
            sampleData(pro_obj)\$sample_id <- sub("\\\\.(bw|bigWig)\$", "", basename(bw_files))

            ht <- generateEnrichedHeatmap(pro_obj, use_raster = FALSE)
            
            pdf("${label}_profile_heatmap.pdf", width=8, height=10)
            print(ht)
            dev.off()

            png_success <- FALSE
            tryCatch({
                png("${label}_profile_heatmap.png", width=1200, height=1400, res=150, type="cairo")
                print(ht)
                dev.off()
                png_success <- TRUE
            }, error = function(e_cairo) {
                png("${label}_profile_heatmap.png", width=1200, height=1400, res=150)
                print(ht)
                dev.off()
                png_success <- TRUE
            })

            if (png_success) {
                img_64 <- base64encode("${label}_profile_heatmap.png")
                cat(paste0(
                    "\\n<div style='text-align: center; padding: 20px;'>\\n",
                    "  <h3>Profile Analysis (", "${label}", ")</h3>\\n",
                    "  <p><strong>Modalità usata: </strong> ", used_mode, " (±2000 bp)</p>\\n",
                    "  <img src='data:image/png;base64,", img_64, "' style='width: 600px; max-width: 100%; height: auto; border: 1px solid #ddd;'>\\n",
                    "</div>"
                ), file="${label}_profileplyr_mqc.html")
            }
        }, error = function(e) {
            cat(paste0(
                "\\n<div style='text-align: center; padding: 20px;'>\\n",
                "  <h3>Profile Analysis: ", "${label}", "</h3>\\n",
                "  <p>Errore durante il plotting: ", e\$message, "</p>\\n",
                "</div>"
            ), file="${label}_profileplyr_mqc.html")
            file.create("${label}_profile_heatmap.png")
            file.create("${label}_profile_heatmap.pdf")
        })
    }

    writeLines(c(
        "\\"${task.process}\\":",
        paste0("    profileplyr: ", packageVersion("profileplyr"))
    ), "versions.yml")
    """
}
