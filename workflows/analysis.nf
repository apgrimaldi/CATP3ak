nextflow.enable.dsl=2

// --- INCLUDE DEI MODULI ---
include { FASTQC }                 from '../modules/local/fastqc.nf'
include { TRIMGALORE }             from '../modules/local/trimgalore.nf'
include { BOWTIE2_BUILD }          from '../modules/local/bowtie2_build.nf'
include { BOWTIE2 }                from '../modules/local/bowtie2.nf'
include { SAMTOOLS_SORT }          from '../modules/local/samtools_sort.nf'
include { SAMTOOLS_STATS }         from '../modules/local/samtools_stats.nf'
include { PICARD_MARKDUPLICATES }  from '../modules/local/picard_markduplicates.nf'
include { FILTERING }              from '../modules/local/filtering.nf'
include { MACS3_ATAC_NARROW }      from '../modules/local/macs3_atac_narrow.nf'
include { MACS3_ATAC_BROAD }       from '../modules/local/macs3_atac_broad.nf'
include { MACS3_CHIP_NARROW }      from '../modules/local/macs3_chip_narrow.nf'
include { MACS3_CHIP_BROAD }       from '../modules/local/macs3_chip_broad.nf'
include { HOMER_ANNOTATEPEAKS }    from '../modules/local/homer_annotate.nf'
include { CALC_FRIP }              from '../modules/local/calc_frip.nf'
include { DEEPTOOLS }              from '../modules/local/deeptools.nf'
include { MULTIQC }                from '../modules/local/multiqc.nf'
include { SAMTOOLS_INDEX }         from '../modules/local/samtools_index.nf'

workflow ATAC_CHIP_PIPELINE {
    take:
    ch_input 

    main:
    ch_versions = Channel.empty()
    ch_multiqc_config = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)

    // --- LOGICA DI ASSEGNAZIONE PARAMETRI GENOMA ---
    def fasta_file     = null
    def gtf_file       = null
    def bowtie2_index  = null
    def blacklist_path = null
    def m_genome       = params.macs_gsize ?: params.genome

    if (params.genomes && params.genomes.containsKey(params.genome)) {
        def gdata = params.genomes[params.genome]
        fasta_file     = params.fasta_file    ?: gdata.fasta
        gtf_file       = params.gtf_file      ?: gdata.gtf
        bowtie2_index  = params.bowtie2_index ?: gdata.bowtie2
        blacklist_path = params.blacklist_file ?: gdata.blacklist
        if (!params.macs_gsize) m_genome = gdata.macs_gsize ?: params.genome
    } else {
        fasta_file     = params.fasta_file
        gtf_file       = params.gtf_file
        bowtie2_index  = params.bowtie2_index
        blacklist_path = params.blacklist_file
    }

    // --- GESTIONE INDICE BOWTIE2 (FIXED) ---
    ch_index_internal = Channel.empty()

    if (bowtie2_index) {
        // Caso 1: L'indice esiste già (path o config)
        ch_index_internal = Channel.fromPath("${bowtie2_index}/*.bt2*").collect()
    } else if (fasta_file) {
        // Caso 2: Abbiamo il file .fna/.fa, costruiamo l'indice
        BOWTIE2_BUILD ( file(fasta_file) )
        // Usiamo .collect() per assicurarci che tutti i file dell'indice vadano insieme
        ch_index_internal = BOWTIE2_BUILD.out.index.collect()
        ch_versions = ch_versions.mix(BOWTIE2_BUILD.out.versions)
    } else {
        error "ERRORE: Impossibile trovare riferimenti per il genoma '${params.genome}'. Fornisci --fasta_file."
    }

    // --- START PIPELINE ---

    // 1. FASTQC
    FASTQC ( ch_input )
    ch_versions = ch_versions.mix(FASTQC.out.versions)

    // 2. TRIMMING
    TRIMGALORE ( ch_input )
    ch_versions = ch_versions.mix(TRIMGALORE.out.versions)

    // 3. ALIGNMENT (FIXED: Passa correttamente il canale dell'indice)
    BOWTIE2 ( TRIMGALORE.out.reads, ch_index_internal )
    ch_versions = ch_versions.mix(BOWTIE2.out.versions)

    // 4. SORTING
    SAMTOOLS_SORT ( BOWTIE2.out.bam )

    // 5. MARK DUPLICATES
    PICARD_MARKDUPLICATES ( SAMTOOLS_SORT.out.bam, [], [] )
    ch_versions = ch_versions.mix(PICARD_MARKDUPLICATES.out.versions)
    
    ch_picard_bams = PICARD_MARKDUPLICATES.out.bam.map { meta, bam, bai -> [ meta, bam, bai ?: [] ] }

    // 6. BLACKLIST FILTERING & INDEXING
    if (blacklist_path) {
        FILTERING ( ch_picard_bams, file(blacklist_path) )
        SAMTOOLS_INDEX ( FILTERING.out.bam )
        ch_final_bams = SAMTOOLS_INDEX.out.bam_bai
        ch_versions = ch_versions.mix(FILTERING.out.versions, SAMTOOLS_INDEX.out.versions)
    } else {
        SAMTOOLS_INDEX ( ch_picard_bams.map { it -> [it[0], it[1]] } )
        ch_final_bams = SAMTOOLS_INDEX.out.bam_bai
    }

    // 7. STATS & DEEPTOOLS
    SAMTOOLS_STATS ( ch_final_bams.map { it -> [ it[0], it[1] ] } )
    ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions)
    DEEPTOOLS ( ch_final_bams )
    ch_versions = ch_versions.mix(DEEPTOOLS.out.versions)

    // 8. PEAK CALLING
    ch_peaks = Channel.empty()
    ch_frip_peaks = Channel.empty()
    ch_macs_logs_mqc = Channel.empty()
    ch_narrow_counts_mqc = Channel.empty()
    ch_broad_counts_mqc  = Channel.empty()

    if (params.protocol == 'atac') {
        ch_macs_input = ch_final_bams.map { it -> [ it[0], it[1] ] }
        MACS3_ATAC_NARROW ( ch_macs_input, m_genome )
        MACS3_ATAC_BROAD  ( ch_macs_input, m_genome )
        
        ch_peaks = MACS3_ATAC_NARROW.out.peaks.mix(MACS3_ATAC_BROAD.out.peaks)
        ch_frip_peaks = MACS3_ATAC_NARROW.out.peaks
        ch_narrow_counts_mqc = MACS3_ATAC_NARROW.out.count_narrow
        ch_broad_counts_mqc  = MACS3_ATAC_BROAD.out.count_broad
        ch_macs_logs_mqc = MACS3_ATAC_NARROW.out.versions.map{ it[1] }.mix(MACS3_ATAC_BROAD.out.versions.map{ it[1] })
    } 
    else if (params.protocol == 'chip') {
        ch_macs3_chip_input = ch_final_bams.map { meta, bam, bai -> [ meta, bam, [] ] } 
        MACS3_CHIP_NARROW ( ch_macs3_chip_input, m_genome )
        MACS3_CHIP_BROAD  ( ch_macs3_chip_input, m_genome )
        
        ch_peaks = MACS3_CHIP_NARROW.out.peaks.mix(MACS3_CHIP_BROAD.out.peaks)
        ch_frip_peaks = MACS3_CHIP_NARROW.out.peaks
        ch_narrow_counts_mqc = MACS3_CHIP_NARROW.out.count_narrow
        ch_broad_counts_mqc  = MACS3_CHIP_BROAD.out.count_broad
        ch_macs_logs_mqc = MACS3_CHIP_NARROW.out.xls.map{ it[1] }.mix(MACS3_CHIP_BROAD.out.xls.map{ it[1] })
    }

    // 9. FRIP & ANNOTATION
    ch_frip_input = ch_final_bams.map { it -> [ it[0], it[1] ] }.join(ch_frip_peaks)
    CALC_FRIP ( ch_frip_input )

    ch_homer_mqc = Channel.empty()
    if (fasta_file && gtf_file) {
        HOMER_ANNOTATEPEAKS ( ch_peaks, file(fasta_file), file(gtf_file) )
        ch_homer_mqc = HOMER_ANNOTATEPEAKS.out.stats.map{ it[1] }.collect().ifEmpty([])
    }

    // 10. MULTIQC
    ch_versions_multiqc = ch_versions.unique().collectFile(name: 'collated_versions.yml')
    ch_all_counts_mqc = ch_narrow_counts_mqc.mix(ch_broad_counts_mqc).map{ it[1] }.collect().ifEmpty([])

    MULTIQC (
        ch_multiqc_config.collect().ifEmpty([]),
        Channel.value("Protocol: ${params.protocol}\nGenome: ${params.genome}\nSize: ${m_genome}").collectFile(name: 'summary.txt'),
        FASTQC.out.zip.map{ it[1] }.collect().ifEmpty([]),
        TRIMGALORE.out.log.map{ it[1] }.collect().ifEmpty([]),
        BOWTIE2.out.log.map{ it[1] }.collect().ifEmpty([]),
        PICARD_MARKDUPLICATES.out.metrics.map{ it[1] }.collect().ifEmpty([]),
        SAMTOOLS_STATS.out.stats.map{ it[1] }.collect().ifEmpty([]),
        DEEPTOOLS.out.fingerprint_txt.map{ it[1] }.collect().ifEmpty([]), 
        ch_macs_logs_mqc.collect().ifEmpty([]), 
        ch_all_counts_mqc,                     
        CALC_FRIP.out.frip.map{ it[1] }.collect().ifEmpty([]), 
        ch_homer_mqc,
        ch_versions_multiqc.collect()                                       
    )
}
