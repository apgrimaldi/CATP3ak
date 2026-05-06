process MACS3_CHIP_BROAD {
    tag "$meta.id"
    label 'process_medium'
    container 'quay.io/biocontainers/macs3:3.0.1--py311h0152c62_3'

    input:
    tuple val(meta), path(ip_bam), path(control_bam)

    output:
    tuple val(meta), path("*.broadPeak") , emit: peaks
    tuple val(meta), path("*.gappedPeak"), emit: gapped_peaks
    path "*.broad_counts.txt"            , emit: count_broad // <--- Output per MultiQC
    path "versions.yml"                  , emit: versions

    script:
    def prefix   = "${meta.id}_broad"
    def format   = meta.single_end ? 'BAM' : 'BAMPE'
    def genome_map = ['hg38': 'hs', 'GRCh38': 'hs', 'hg19': 'hs', 'mm10': 'mm', 'mm9': 'mm', 'GRCm38': 'mm', 'dm6': 'dm', 'ce11': 'ce']
    def m_genome = genome_map[params.genome] ?: params.genome

    """
    macs3 callpeak \\
        -t $ip_bam \\
        -c $control_bam \\
        -f $format \\
        -g $m_genome \\
        -n $prefix \\
        --broad \\
        --broad-cutoff 0.1

    # Genera il file di conteggio specifico per i BROAD
    count=\$(wc -l < ${prefix}.broadPeak)
    echo -e "${meta.id}\t\$count" > ${meta.id}.broad_counts.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        macs3: \$(macs3 --version | sed 's/macs3 //g')
    END_VERSIONS
    """
}
