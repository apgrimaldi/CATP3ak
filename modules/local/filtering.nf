process FILTERING {
    tag "$meta.id"
    label 'process_medium'
    
    container 'quay.io/biocontainers/bedtools:2.30.0--hc088bd4_0'

    publishDir "${params.outdir}/04_filtered", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)
    path  blacklist

    output:
    tuple val(meta), path("*.filtered.bam"), emit: bam
    path "versions.yml"                    , emit: versions

    script:
    def prefix = "${meta.id}"
    """
    # 1. Gestione Blacklist: 
    # Se il file è .gz lo decomprimiamo, altrimenti creiamo un link simbolico
    if [[ "$blacklist" == *.gz ]]; then
        gunzip -c "$blacklist" > actual_blacklist.bed
    else
        # Usiamo ln -s per non duplicare i dati se la blacklist è già un .bed
        ln -s $blacklist actual_blacklist.bed
    fi

    # 2. Rimozione Blacklist con bedtools
    # -v: 
    # -abam: 
    bedtools intersect \\
        -v \\
        -abam $bam \\
        -b actual_blacklist.bed \\
        > ${prefix}.filtered.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed 's/bedtools v//')
    END_VERSIONS
    """
}
