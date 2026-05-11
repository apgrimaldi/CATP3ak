process DEEPTOOLS {
    tag "$meta.id"
    label 'process_high'
    container 'quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0'
    
    publishDir "${params.outdir}/07_advanced_qc", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.bigWig")                           , emit: bw
    tuple val(meta), path("*.plotFingerprint.pdf")              , emit: fingerprint_pdf
    tuple val(meta), path("*.plotFingerprint.raw.txt")          , emit: fingerprint_txt 
    tuple val(meta), path("*.plotFingerprint.qcmetrics.txt")    , emit: fingerprint_metrics
    path "versions.yml"                                         , emit: versions

    script:
    def prefix = "${meta.id}"
    // Gestione reads single-end come fa nf-core
    def extend = (meta.single_end && params.fragment_size > 0) ? "--extendReads ${params.fragment_size}" : ''
    """
    # 1. Genera BigWig
    bamCoverage \\
        -b $bam \\
        -o ${prefix}.bigWig \\
        --binSize 10 \\
        --normalizeUsing CPM \\
        --numberOfProcessors $task.cpus

    # 2. Fingerprint (Parametri nf-core)
    plotFingerprint \\
        --bamfiles $bam \\
        --plotFile ${prefix}.plotFingerprint.pdf \\
        --outRawCounts ${prefix}.plotFingerprint.raw.txt \\
        --outQualityMetrics ${prefix}.plotFingerprint.qcmetrics.txt \\
        --numberOfProcessors $task.cpus \\
        $extend \\
        --skipZeros

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: \$(plotFingerprint --version | sed -e "s/plotFingerprint //g")
    END_VERSIONS
    """
}
