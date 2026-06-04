process BOWTIE2 {
    tag "$meta.id"
    label 'process_high'
    container 'quay.io/biocontainers/mulled-v2-c742dccc9d8fabfcff2af0d8d6799dbc711366cf:0b0dea2b5dffed0cff6fb77b4377a5940cd4319a-0'

    input:
    tuple val(meta), path(reads)
    path index

    output:
    tuple val(meta), path("*.raw.bam"), emit: bam
    tuple val(meta), path("*.log")    , emit: log
    path "versions.yml"               , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def rg_args = "--rg-id ${prefix} --rg SM:${prefix} --rg PL:ILLUMINA --rg LB:lib1"
    
    def input_reads = meta.single_end ? "-U ${reads}" : "-1 ${reads[0]} -2 ${reads[1]}"
    
    def extra_args = params.protocol == 'atac' ? "--no-mixed --no-discordant" : ""

    """
    INDEX_BASE=\$(find -L . -name "*.1.bt2*" | sed 's/\\.1\\.bt2.*//' | head -n 1)

    if [ -z "\$INDEX_BASE" ]; then
        echo "ERRORE: Indice Bowtie2 non trovato nella cartella di lavoro."
        ls -la
        exit 1
    fi

    bowtie2 \\
        -x "\$INDEX_BASE" \\
        $input_reads \\
        -p $task.cpus \\
        $rg_args \\
        --very-sensitive \\
        $extra_args \\
        -X 2000 \\
        2> ${prefix}.bowtie2.log \\
        | samtools view -@ $task.cpus -b -o ${prefix}.raw.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(echo \$(bowtie2 --version 2>&1) | sed 's/^.*bowtie2-align-s version //; s/ .*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
