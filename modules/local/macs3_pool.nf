process MACS3_POOL {
tag "$meta.id"
label 'process_high'
container 'quay.io/biocontainers/macs3:3.0.0a7--py310h5624777_1'

input:
tuple val(meta), path(ip_bams), path(control_bams)
val macs_gsize

output:
tuple val(meta), path("*_peaks.narrowPeak"), emit: peaks
tuple val(meta), path("*_peaks.xls")       , emit: xls
path "*_peaks.xls"                         , emit: count_narrow
path "versions.yml"                        , emit: versions

script:
// Qui definiamo i parametri avanzati che hai richiesto.
// In futuro, puoi anche sovrascriverli da nextflow.config usando ext.args
def args = task.ext.args ?: '--nomodel --shift -100 --extsize 200 -B'

"""
# Nextflow trasforma automaticamente le liste "ip_bams" e "control_bams" in 
# stringhe separate da spazio (es: "file1.bam file2.bam"), che e' il formato 
# esatto richiesto da MACS3 per fare il pooling!

macs3 callpeak \\
    -t ${ip_bams} \\
    -c ${control_bams} \\
    -f BAM \\
    -g ${macs_gsize} \\
    -n ${meta.id} \\${args}

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    macs3: \$(macs3 --version | sed -e "s/macs3 //g")
END_VERSIONS
"""


}
