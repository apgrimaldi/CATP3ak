process DEEPTOOLS {
    tag "$meta.id"
    label 'process_high'
    container 'quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0'
    
    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.display.bw")           , emit: bw_display
    tuple val(meta), path("*.l6n.bw")               , emit: bw_lanceotron
    tuple val(meta), path("*.plotFingerprint.pdf")  , emit: fingerprint_pdf
    tuple val(meta), path("*.plotFingerprint.raw.txt"), emit: fingerprint_txt 
    tuple val(meta), path("*.plotFingerprint.qcmetrics.txt"), emit: fingerprint_metrics
    path "versions.yml"                             , emit: versions

    script:
    def prefix = "${meta.id}"
    def extend = (meta.single_end && params.fragment_size > 0) ? "--extendReads ${params.fragment_size}" : (params.single_end ? "--extendReads" : "")
    
    """
   
    bamCoverage \\
        --bam $bam \\
        --outFileName ${prefix}.display.bw \\
        --binSize 10 \\
        --normalizeUsing CPM \\
        --numberOfProcessors $task.cpus \\
        $extend

  
    bamCoverage \\
        --bam $bam \\
        --outFileName ${prefix}.l6n.bw \\
        --binSize 1 \\
        --normalizeUsing RPKM \\
        --numberOfProcessors $task.cpus \\
        $extend

 
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
