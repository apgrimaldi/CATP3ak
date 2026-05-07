workflow {
    
    // 1. Controllo validità
    if (!params.input) { error "Errore: Specifica --input" }
    
    log.info """
    ===========================================
    P I P E L I N E   A T A C / C H I P
    ===========================================
    Protocollo : ${params.protocol?.toUpperCase()}
    Genoma     : ${params.genome}
    Input      : ${params.input}
    ===========================================
    """

    // 2. Lettura Samplesheet
    ch_input = Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header:true, sep:',')
        .map { row -> create_fastq_channel(row) }

    // 3. Lancio della Pipeline
    // Abbiamo spostato tutta la logica dell'indice dentro ATAC_CHIP_PIPELINE
    // Passiamo solo il canale degli input (ch_input)
    ATAC_CHIP_PIPELINE ( ch_input )
    
    workflow.onComplete {
        log.info "Pipeline completata con successo!"
    }
}
