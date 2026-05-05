nextflow.enable.dsl=2

include { ATAC_CHIP_PIPELINE } from './workflows/analysis'

def create_fastq_channel(LinkedHashMap row) {
    def meta = [:]
    meta.id         = row.sample.trim()
    meta.antibody   = (row.antibody && row.antibody.trim() != "") ? row.antibody.trim() : 'none'
    meta.control    = (row.control && row.control.trim() != "") ? row.control.trim() : 'none'
    
    // Verifica Single-End o Paired-End
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

    // 3. Gestione Indice Bowtie2 (Corretta per leggere i file .bt2)
    def index_path = params.bowtie2_index ?: params.genomes[ params.genome ]?.bowtie2 ?: null

    if (!index_path) {
        error "Errore: Indice Bowtie2 non trovato per '${params.genome}'."
    }

    // Carichiamo tutti i file dell'indice (.bt2 o .bt2l)
    ch_index = Channel
        .fromPath("${index_path}/*.bt2*", checkIfExists: true)
        .collect()

    // 4. Lancio della Pipeline
    // IMPORTANTE: Assicurati che 'analysis.nf' accetti questi due argomenti in 'take'
    ATAC_CHIP_PIPELINE ( ch_input, ch_index )
    
    workflow.onComplete {
        log.info "Pipeline completata con successo!"
    }
}
