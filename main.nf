nextflow.enable.dsl=2

include { ATAC_CHIP_PIPELINE } from './workflows/analysis' 

def create_fastq_channel(LinkedHashMap row, Set known_controls) {
    def meta = [:]
    meta.id         = row.sample.trim()
    meta.antibody   = (row.antibody && row.antibody.trim() != "") ? row.antibody.trim() : 'none'
    meta.control    = (row.control && row.control.trim() != "") ? row.control.trim() : 'none'
    
    if (params.protocol == 'atac') {
        meta.is_control = false
    } else {
        meta.is_control = known_controls.contains(meta.id) || meta.antibody.toLowerCase() == 'igg' || row.is_control == 'true'
    }

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
    if (!params.input) { error "Error: Please specify --input samplesheet.csv" }
    
    def known_controls = [] as Set
    file(params.input).splitCsv(header:true, sep:',').each { row ->
        if (row.control && row.control.trim() != "") {
            known_controls.add(row.control.trim())
        }
    }

    log.info """
    ===========================================
    A T A C / C H I P   P I P E L I N E
    ===========================================
    Protocol      : ${params.protocol?.toUpperCase()}
    Genome        : ${params.genome}
    Input         : ${params.input}
    Output        : ${params.outdir}
    ===========================================
    """

    ch_input = Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header:true, sep:',')
        .map { row -> create_fastq_channel(row, known_controls) }
    
    ch_input.view { meta, reads -> 
        "LOG: ID: ${meta.id} | Group: ${meta.antibody} | Control: ${meta.is_control}" 
    }

    ATAC_CHIP_PIPELINE ( ch_input )
    
    workflow.onComplete {
        log.info "Pipeline completed successfully!"
    }
}
