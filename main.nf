nextflow.enable.dsl=2

include { CATP3ak } from './workflows/analysis' 

// --- HELP MESSAGE FUNCTION ---
def helpMessage() {
    log.info"""
    =======================================================================
                   C A T P 3 A K   P I P E L I N E
    =======================================================================
    A Nextflow DSL2 pipeline for comprehensive ATAC-seq and ChIP-seq 
    data analysis.

    Usage:
    nextflow run main.nf --input samplesheet.csv --protocol chip --genome hg38

    -----------------------------------------------------------------------
    MANDATORY ARGUMENTS:
    -----------------------------------------------------------------------
      --input                  Path to comma-separated file containing sample 
                               information (CSV). Required columns: 
                               sample, fastq_1, fastq_2, antibody, control, group.
      --protocol               Specify the experimental protocol. 
                               Valid options: 'chip' or 'atac'. (Default: chip)
      --genome                 Reference genome abbreviation (e.g., 'hg38', 'mm10').

    -----------------------------------------------------------------------
    GENERAL PARAMETERS:
    -----------------------------------------------------------------------
      --outdir                 The output directory where results will be saved. 
                               (Default: 'results')
      --single_end             Set to 'true' if the input reads are single-end. 
                               (Default: false)
      --fragment_size          Estimated fragment size used for extending reads 
                               in DeepTools. (Default: 200)
      --lanceotron_threshold   Confidence threshold for filtering Lanceotron peaks. 
                               (Default: 0.5)

    -----------------------------------------------------------------------
    CUSTOM REFERENCE & GENOME OPTIONS (Overrides --genome defaults):
    -----------------------------------------------------------------------
      --reference_file         Path to the FASTA reference genome.
      --gtf_file               Path to the GTF annotation file (used by HOMER).
      --bowtie2_index          Path to the pre-built Bowtie2 index directory.
      --blacklist              Path to the BED file containing blacklisted regions.
      --chrom_sizes            Path to the chromosome sizes file (Required for Omnipeak).
      --macs_gsize             Effective genome size for MACS3 (e.g., 'hs', 'mm').

    -----------------------------------------------------------------------
    SKIP OPTIONS (Disable specific analysis steps):
    -----------------------------------------------------------------------
      --skip_lanceotron        Skip peak calling using Lanceotron.
      --skip_macs3             Skip peak calling using MACS3.
      --skip_omnipeak          Skip peak calling using Omnipeak.
      --skip_frip              Skip the calculation of FRiP (Fraction of Reads in Peaks).
      --skip_homer             Skip peak annotation using HOMER.
      --skip_diffbind          Skip differential binding analysis using DiffBind.
      --skip_profileplyr       Skip signal profile generation using ProfilePlyr.

    -----------------------------------------------------------------------
    SYSTEM OPTIONS:
    -----------------------------------------------------------------------
      --help                   Display this help message and exit.
    =======================================================================
    """.stripIndent()
}

// --- CHECK FOR HELP FLAG ---
if (params.help) {
    helpMessage()
    exit 0
}

// --- CHANNEL PREPARATION FUNCTION ---
def create_fastq_channel(LinkedHashMap row, Set known_controls) {
    def meta = [:]
    meta.id         = row.sample.trim()
    meta.antibody   = (row.antibody && row.antibody.trim() != "") ? row.antibody.trim() : 'none'
    meta.control    = (row.control && row.control.trim() != "") ? row.control.trim() : 'none'
    meta.group      = (row.group && row.group.trim() != "") ? row.group.trim() : 'Baseline'
    
    // Automatically flag controls based on protocol and samplesheet data
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

// --- MAIN WORKFLOW EXECUTION ---
workflow {
    if (!params.input) { 
        error "Error: Please specify an input samplesheet using --input samplesheet.csv or run 'nextflow run apgrimaldi/CATP3ak - latest --help' for available options." 
    }
    
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
    
    // Print input information
    ch_input.view { meta, reads -> 
        "LOG: ID: ${meta.id} | Antibody: ${meta.antibody} | Group: ${meta.group} | Control: ${meta.is_control}" 
    }

    CATP3ak ( ch_input )
    
    workflow.onComplete {
        log.info "CATP3ak Pipeline completed successfully!"
    }
}
