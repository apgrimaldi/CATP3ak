nextflow.enable.dsl=2

// --- INCLUDE DEI MODULI ---
// (Mantieni i tuoi include invariati...)
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
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_FINAL } from '../modules/local/samtools_index.nf'

// NUOVO INCLUDE: DiffBind
include { DIFFBIND }               from '../modules/local/diffbind.nf'

workflow ATAC_CHIP_PIPELINE {
    take:
    ch_input 

    main:
    ch_versions = Channel.empty()
    ch_multiqc_config = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)

    // --- 1. LOGICA GENOMA (Invariata) ---
    // ... [Codice genoma omesso per brevità] ...

    // --- 2. INDICE BOWTIE2 (Invariata) ---
    // ... [Codice bowtie2 build omesso per brevità] ...

    // --- 3. CORE PROCESSING ---
    FASTQC ( ch_input )
    TRIMGALORE ( ch_input )
    ch_versions = ch_versions.mix(FASTQC.out.versions, TRIMGALORE.out.versions)

    BOWTIE2 ( TRIMGALORE.out.reads, ch_index_internal )
    ch_versions = ch_versions.mix(BOWTIE2.out.versions)

    SAMTOOLS_SORT ( BOWTIE2.out.bam )

    def ch_fasta_ref = fasta_file ? file(fasta_file) : []
    PICARD_MARKDUPLICATES ( SAMTOOLS_SORT.out.bam, ch_fasta_ref, [] )
    ch_versions = ch_versions.mix(PICARD_MARKDUPLICATES.out.versions)

    SAMTOOLS_INDEX ( PICARD_MARKDUPLICATES.out.bam )
    
    if (blacklist_path) {
        FILTERING ( SAMTOOLS_INDEX.out.bam_bai, file(blacklist_path) )
        SAMTOOLS_INDEX_FINAL ( FILTERING.out.bam )
        ch_final_bams = SAMTOOLS_INDEX_FINAL.out.bam_bai
    } else {
        ch_final_bams = SAMTOOLS_INDEX.out.bam_bai
    }

    // --- 4. METRICHE & PEAKS ---
    SAMTOOLS_STATS ( ch_final_bams.map { meta, bam, bai -> [ meta, bam ] } )
    DEEPTOOLS ( ch_final_bams )

    ch_macs_input = ch_final_bams.map { meta, bam, bai -> [ meta, bam ] }
    
    ch_peaks = Channel.empty()
    ch_frip_peaks = Channel.empty()
    // ... [Inizializzazione canali macs invariata] ...

    if (params.protocol == 'atac') {
        MACS3_ATAC_NARROW ( ch_macs_input, m_genome )
        MACS3_ATAC_BROAD ( ch_macs_input, m_genome )
        ch_peaks = MACS3_ATAC_NARROW.out.peaks.mix(MACS3_ATAC_BROAD.out.peaks)
        ch_frip_peaks = MACS3_ATAC_NARROW.out.peaks
        // ... [mapping canali macs atac] ...
    } else {
        ch_macs_chip_input = ch_final_bams.map { meta, bam, bai -> [ meta, bam, [] ] }
        MACS3_CHIP_NARROW ( ch_macs_chip_input, m_genome )
        MACS3_CHIP_BROAD ( ch_macs_chip_input, m_genome )
        ch_peaks = MACS3_CHIP_NARROW.out.peaks.mix(MACS3_CHIP_BROAD.out.peaks)
        ch_frip_peaks = MACS3_CHIP_NARROW.out.peaks
        // ... [mapping canali macs chip] ...
    }

    // --- 5. ANNOTAZIONE E FRIP ---
    ch_frip_input = ch_final_bams.map { meta, bam, bai -> [ meta, bam ] }.join(ch_frip_peaks)
    CALC_FRIP ( ch_frip_input )

    ch_homer_mqc = Channel.empty()
    if (fasta_file && gtf_file) {
        HOMER_ANNOTATEPEAKS ( ch_peaks, file(fasta_file), file(gtf_file) )
        ch_homer_mqc = HOMER_ANNOTATEPEAKS.out.stats_mqc.map{ it[1] }.collect().ifEmpty([])
    }

    // --- 6. DIFFBIND (NUOVA SEZIONE) ---
    ch_diffbind_mqc = Channel.empty()
    if (params.samplesheet_diffbind) {
        // Raccogliamo tutti i file necessari per DiffBind
        def ch_bams_db  = ch_final_bams.map { it[1] }.collect()
        def ch_bais_db  = ch_final_bams.map { it[2] }.collect()
        def ch_peaks_db = ch_peaks.map { it[1] }.collect()

        DIFFBIND (
            file(params.samplesheet_diffbind),
            ch_bams_db,
            ch_bais_db,
            ch_peaks_db
        )
        ch_diffbind_mqc = DIFFBIND.out.mqc_html.collect().ifEmpty([])
        ch_versions = ch_versions.mix(DIFFBIND.out.versions)
    }

    // --- 7. MULTIQC ---
    ch_versions_multiqc = ch_versions.unique().collectFile(name: 'collated_versions.yml')
    ch_all_counts_mqc = ch_narrow_counts_mqc.mix(ch_broad_counts_mqc).map{ it[1] }.collect().ifEmpty([])

    MULTIQC (
        ch_multiqc_config.collect().ifEmpty([]),                                      
        Channel.value("Protocol: ${params.protocol}\nGenome: ${params.genome}").collectFile(name: 'summary.txt'), 
        FASTQC.out.zip.map{ it[1] }.collect().ifEmpty([]),                            
        TRIMGALORE.out.log.map{ it[1] }.collect().ifEmpty([]),                        
        BOWTIE2.out.log.map{ it[1] }.collect().ifEmpty([]),                           
        PICARD_MARKDUPLICATES.out.metrics.map{ it[1] }.collect().ifEmpty([]),         
        SAMTOOLS_STATS.out.stats.map{ it[1] }.collect().ifEmpty([]),                  
        DEEPTOOLS.out.fingerprint_txt.map{ it[1] }.mix(DEEPTOOLS.out.fingerprint_metrics.map{ it[1] }).collect().ifEmpty([]),
        ch_macs_logs_mqc.collect().ifEmpty([]),                                       
        ch_all_counts_mqc,                                                            
        CALC_FRIP.out.frip.map{ it[1] }.collect().ifEmpty([]),                        
        ch_homer_mqc,                                                                 
        ch_diffbind_mqc, // AGGIUNTO: Canale HTML di DiffBind
        ch_versions_multiqc.collect()                                                 
    )
}
