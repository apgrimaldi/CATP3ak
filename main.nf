nextflow.enable.dsl=2

// Importa il workflow con il SUO NUOVO NOME dal file analysis.nf
include { CATP3ak } from './workflows/analysis' 

def create_fastq_channel(LinkedHashMap row, Set known_controls) {
    def meta = [:]
    meta.id         = row.sample.trim()
    meta.antibody   = (row.antibody && row.antibody.trim() != "") ? row.antibody.trim() : 'none'
    meta.control    = (row.control && row.control.trim() != "") ? row.control.trim() : 'none'
    meta.group      = (row.group && row.group.trim() != "") ? row.group.trim() : 'Baseline'
    
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
         C A T P 3 A K   P I P E L I N E
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
        "LOG: ID: ${meta.id} | Antibody: ${meta.antibody} | Group: ${meta.group} | Control: ${meta.is_control}" 
    }

    // Richiama il workflow usando il nuovo acronimo
    CATP3ak ( ch_input )
    
    workflow.onComplete {
        log.info "CATP3ak Pipeline completed successfully!"
    }
}
