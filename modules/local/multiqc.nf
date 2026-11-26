process MULTIQC {
    label 'process_medium'
    container 'quay.io/biocontainers/multiqc:1.34--pyhdfd78af_0'

    input:
    path multiqc_config
    path workflow_summary
    path fastqc
    path trimgalore
    path alignment
    path picard
    path samtools
    path deeptools
    path macs3
    path counts
    path frip
    path homer
    path diffbind
    path seqcode
    path deeptools_heatmaps
    path lanceotron
    path omnipeak
    path versions

    output:
    path "*multiqc_report.html", emit: report
    path "*_data"               , emit: data
    path "versions.yml"         , emit: versions

    script:
    def args = task.ext.args ?: ''
    def config_opt = (multiqc_config && multiqc_config.name != 'empty_config' && multiqc_config.size() > 0) ? "--config $multiqc_config" : ''
    def report_title = params.protocol == 'atac' ? "ATAC-seq Analysis Report" : "ChIP-seq Analysis Report"

    """
    python3 - << 'EOF'
    import base64
    import glob
    import html

    # 1) Sezione dedicata alle heatmap di deepTools (*.png)
    pngs = sorted(glob.glob("*.png"))
    if pngs:
        with open("deeptools_heatmaps_mqc.html", "w", encoding="utf-8") as out:
            out.write("<!--\\n")
            out.write("id: deeptools_tss_heatmaps\\n")
            out.write("section_name: 'DeepTools: TSS Heatmaps'\\n")
            out.write("description: 'Heatmap del segnale attorno ai TSS generate con deepTools plotHeatmap.'\\n")
            out.write("-->\\n")
            for p in pngs:
                with open(p, "rb") as f:
                    b64 = base64.b64encode(f.read()).decode("ascii")
                out.write(
                    f'<div style="margin-bottom: 24px;">'
                    f'<h4>{html.escape(p)}</h4>'
                    f'<img src="data:image/png;base64,{b64}" style="max-width: 100%; height: auto;" />'
                    f'</div>\\n'
                )

    # 2) Sezione dedicata a SeqCode genomeDistribution
    seqcode_files = sorted(glob.glob("*_seqcode*"))
    if seqcode_files:
        with open("seqcode_dist_mqc.html", "w", encoding="utf-8") as out:
            out.write("<!--\\n")
            out.write("id: seqcode_genome_distribution\\n")
            out.write("section_name: 'SeqCode: Genome Distribution'\\n")
            out.write("description: 'Distribuzione genomica dei picchi calcolata con SeqCode genomeDistribution.'\\n")
            out.write("-->\\n")
            for sf in seqcode_files:
                if sf.endswith((".png", ".jpg", ".jpeg", ".pdf", ".svg")):
                    continue
                with open(sf, "r", encoding="utf-8", errors="replace") as f:
                    content = html.escape(f.read())
                out.write(
                    f'<div style="margin-bottom: 20px;">'
                    f'<h4>{html.escape(sf)}</h4>'
                    f'<pre style="max-height: 360px; overflow: auto;">{content}</pre>'
                    f'</div>\\n'
                )
    EOF

    multiqc \\
        -f \\
        $args \\
        $config_opt \\
        --title "$report_title" \\
        .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$(multiqc --version | sed 's/multiqc, version //g')
    END_VERSIONS
    """
}
