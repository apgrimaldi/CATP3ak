process TRIMGALORE {
    tag "${meta.id}"
    label 'process_high' 

    container 'quay.io/biocontainers/trim-galore:0.6.11--hdfd78af_0'

    input:
    tuple val(meta), path(reads)

    output:

    tuple val(meta), path("*.fq.gz")       , emit: reads
    tuple val(meta), path("*_report.txt")  , emit: log
    path "versions.yml"                    , emit: versions

    script:
    
    def cores = task.cpus ? Math.max(Math.floor(task.cpus / 2) as int, 1) : 2

    if (meta.single_end) {
        
        """
        trim_galore \\
            --cores $cores \\
            --gzip \\
            $reads

        cat <<EOF > versions.yml
        "${task.process}":
            trimgalore: \$(echo \$(trim_galore --version 2>&1) | sed 's/^.*version //; s/ .*\$//')
            cutadapt: \$(cutadapt --version | head -n 1)
        EOF
        """
    } else {
        
        """
        trim_galore \\
            --cores $cores \\
            --paired \\
            --gzip \\
            ${reads[0]} \\
            ${reads[1]}

        cat <<EOF > versions.yml
        "${task.process}":
            trimgalore: \$(echo \$(trim_galore --version 2>&1) | sed 's/^.*version //; s/ .*\$//')
            cutadapt: \$(cutadapt --version | head -n 1)
        EOF
        """
    }
}
