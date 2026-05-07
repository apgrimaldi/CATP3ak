nextflow.enable.dsl=2

// --- IMPORTAZIONE WORKFLOW ---
// Se analysis.nf è dentro la cartella 'workflows', usa questo:
include { ATAC_CHIP_PIPELINE } from './workflows/analysis' 
// Se invece analysis.nf è nella stessa cartella di main.nf, usa: include { ATAC_CHIP_PIPELINE } from './analysis'

def create_fastq_channel(LinkedHashMap row) {
    def meta = [:]
    meta.id         = row.sample.trim()
    meta.antibody   = (row.antibody && row.antibody.trim() != "") ? row.antibody.trim() : 'none'
    meta.control    = (row.control && row.control.trim() != "") ? row.control.trim() : 'none'
    meta.single_end = (row.fastq_2 == null || row.fastq_2.trim() == "") ? true : false

    def fastq_1 = file(row.fastq_1, checkIfExists: true)
    def fastqs = [ fastq_1 ]
    
    if (!meta.single_end) {
        def fastq_2 = file(row.fastq_2, checkIfExists: true)
        fastqs << fastq_2
    }
    return [ meta, fastqs ]
}

workflow {
    // 1. Controllo validità input
    if (!params.input) { error "Errore: Specifica --input samplesheet.csv" }
    
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
    // Passiamo il canale al workflow definito in analysis.nf
    ATAC_CHIP_PIPELINE ( ch_input )
    
    workflow.onComplete {
        log.info "Pipeline completata con successo!"
    }
}
