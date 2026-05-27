process PROFILEPLYR {
    tag "${label}"
    label 'process_high'
    container 'quay.io/biocontainers/bioconductor-profileplyr:1.22.0--r44hdfd78af_0'

    input:
    path peaks
    path bigwigs
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

    # 1. Identificazione corretta dei file passati da Nextflow
    peak_files <- list.files(pattern = "\\\\.(bed|narrowPeak|broadPeak)\$")
    bw_files <- list.files(pattern = "\\\\.(bw|bigWig)\$")

    print("File di picchi trovati:")
    print(peak_files)
    print("File BigWig trovati:")
    print(bw_files)

    # Funzione di sicurezza per unire i picchi se sono multipli
    import_and_merge_peaks <- function(files) {
        gr_list <- lapply(files, function(f) {
            tryCatch({
                gr <- rtracklayer::import(f)
                # Forza lo score a numerico per evitare conflitti nell'unione
                if (!is.null(mcols(gr)\$score)) mcols(gr)\$score <- as.numeric(mcols(gr)\$score)
                return(gr)
            }, error = function(e) NULL)
        })
        gr_list <- gr_list[!sapply(gr_list, is.null)]
        if (length(gr_list) == 0) return(NULL)
        
        # Unisce le regioni riducendo quelle sovrapposte
        merged <- GenomicRanges::reduce(do.call(c, gr_list))
        return(merged)
    }

    peaks_gr <- import_and_merge_peaks(peak_files)

    # Controllo di sicurezza: se non ci sono picchi utili, non crashare!
    if (is.null(peaks_gr) || length(peaks_gr) == 0) {
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
        
        # 2. Creazione corretta dell'oggetto Profileplyr per BigWig
        # Esportiamo i picchi uniti su file BED temporaneo (BamBigwig_to_chipProfile lo richiede fisicamente)
        rtracklayer::export(peaks_gr, "merged_testRanges.bed")

        # Usiamo la VERA funzione del pacchetto per leggere i file BigWig
        pro_chip <- BamBigwig_to_chipProfile(
            signalFiles = bw_files,
            testRanges = "merged_testRanges.bed",
            format = "bigwig",
            style = "point",
            bin_size = 50,
            distanceAround = 2000
        )

        # Convertiamo l'output nell'oggetto profileplyr standard
        pro_obj <- as_profileplyr(pro_chip)

        # Pulizia dei nomi dei campioni per la legenda
        sampleData(pro_obj)\$sample_id <- sub("\\\\.(bw|bigWig)\$", "", basename(bw_files))

        # 3. Generazione e salvataggio Heatmap grafica
        tryCatch({
            png("${label}_profile_heatmap.png", width=1200, height=1400, res=150, type="cairo")
            ht <- generateEnrichedHeatmap(pro_obj)
            print(ht)
            dev.off()

            pdf("${label}_profile_heatmap.pdf", width=8, height=10)
            print(ht)
            dev.off()

            # 4. Preparazione HTML codificato in Base64 per MultiQC
            img_64 <- base64encode("${label}_profile_heatmap.png")
            cat(paste0(
                "\\n",
                "<div style='text-align: center; padding: 20px;'>\\n",
                "  <h3>Profile Analysis (", "${label}", ")</h3>\\n",
                "  <img src='data:image/png;base64,", img_64, "' style='width: 600px; max-width: 100%; height: auto; border: 1px solid #ddd;'>\\n",
                "</div>"
            ), file="${label}_profileplyr_mqc.html")

        }, error = function(e) {
            cat(paste0(
                "\\n",
                "<div style='text-align: center; padding: 20px;'>\\n",
                "  <h3>Profile Analysis: ", "${label}", "</h3>\\n",
                "  <p>Errore durante il plotting: ", e\$message, "</p>\\n",
                "</div>"
            ), file="${label}_profileplyr_mqc.html")
            file.create("${label}_profile_heatmap.png")
            file.create("${label}_profile_heatmap.pdf")
        })
    }

    # Scrittura versioni software
    writeLines(c(
        "\\"${task.process}\\":",
        paste0("    profileplyr: ", packageVersion("profileplyr")),
        paste0("    rtracklayer: ", packageVersion("rtracklayer"))
    ), "versions.yml")
    """
}
