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
    path seqcode, stageAs: 'seqcode_in/*'
    path deeptools_heatmaps, stageAs: 'heatmaps_in/*'
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
    import os

    def build_carousel(carousel_id, title_text, items):
        # items: list of dicts with 'label' and 'inner_html'
        slides_html = []
        for idx, it in enumerate(items):
            disp = "block" if idx == 0 else "none"
            slides_html.append(
                f'<div class="{carousel_id}_slide" style="display:{disp}; text-align:center;">'
                f'<div style="font-weight:600; margin-bottom:8px;">{html.escape(it["label"])} ({idx + 1}/{len(items)})</div>'
                f'{it["inner_html"]}'
                f'</div>'
            )

        controls = ""
        script = ""
        if len(items) > 1:
            controls = (
                f'<div style="display:flex; justify-content:center; align-items:center; gap:16px; margin-bottom:12px;">'
                f'<button type="button" onclick="move_{carousel_id}(-1)" style="padding:6px 14px; cursor:pointer; font-size:16px;">&#9664; Prev</button>'
                f'<span id="{carousel_id}_counter" style="font-weight:600;">1 / {len(items)}</span>'
                f'<button type="button" onclick="move_{carousel_id}(1)" style="padding:6px 14px; cursor:pointer; font-size:16px;">Next &#9654;</button>'
                f'</div>'
            )
            script = (
                f'<script>'
                f'(function() {{'
                f'  let idx = 0;'
                f'  const slides = document.getElementsByClassName("{carousel_id}_slide");'
                f'  const counter = document.getElementById("{carousel_id}_counter");'
                f'  window.move_{carousel_id} = function(step) {{'
                f'    if (!slides.length) return;'
                f'    slides[idx].style.display = "none";'
                f'    idx = (idx + step + slides.length) % slides.length;'
                f'    slides[idx].style.display = "block";'
                f'    if (counter) counter.textContent = (idx + 1) + " / " + slides.length;'
                f'  }};'
                f'}})();'
                f'</script>'
            )

        return (
            f'<div style="border:1px solid #ddd; border-radius:8px; padding:16px; margin:12px 0;">'
            f'{controls}'
            f'{"".join(slides_html)}'
            f'</div>'
            f'{script}'
        )

    # 1) DeepTools TSS Heatmaps (*.png)
    pngs = sorted(glob.glob("heatmaps_in/*.png"))
    if pngs:
        dt_items = []
        for p in pngs:
            fname = os.path.basename(p)
            with open(p, "rb") as f:
                b64 = base64.b64encode(f.read()).decode("ascii")
            dt_items.append({
                "label": fname,
                "inner_html": f'<img src="data:image/png;base64,{b64}" style="max-width:100%; max-height:720px; height:auto;" />'
            })
        with open("deeptools_heatmaps_mqc.html", "w", encoding="utf-8") as out:
            out.write("<!--\\n")
            out.write("id: deeptools_tss_heatmaps\\n")
            out.write("section_name: 'DeepTools: TSS Heatmaps'\\n")
            out.write("description: 'Heatmap del segnale attorno ai TSS (scorribili con le frecce).'\\n")
            out.write("-->\\n")
            out.write(build_carousel("dt_hm_carousel", "DeepTools TSS Heatmaps", dt_items))

    # 2) SeqCode Genome Distribution (immagini se presenti, altrimenti testo/tabelle)
    sc_all = sorted([f for f in glob.glob("seqcode_in/*") if os.path.isfile(f)])
    sc_imgs = [f for f in sc_all if f.lower().endswith((".png", ".jpg", ".jpeg", ".webp"))]
    sc_texts = [f for f in sc_all if not f.lower().endswith((".png", ".jpg", ".jpeg", ".webp", ".pdf", ".svg"))]

    sc_items = []
    for img_path in sc_imgs:
        fname = os.path.basename(img_path)
        ext = os.path.splitext(fname)[1].lower().lstrip(".")
        mime = "jpeg" if ext == "jpg" else ext
        with open(img_path, "rb") as f:
            b64 = base64.b64encode(f.read()).decode("ascii")
        sc_items.append({
            "label": fname,
            "inner_html": f'<img src="data:image/{mime};base64,{b64}" style="max-width:100%; max-height:720px; height:auto;" />'
        })

    for txt_path in sc_texts:
        fname = os.path.basename(txt_path)
        with open(txt_path, "r", encoding="utf-8", errors="replace") as f:
            content = html.escape(f.read())
        sc_items.append({
            "label": fname,
            "inner_html": f'<pre style="text-align:left; max-height:420px; overflow:auto; background:#f8f9fa; padding:12px;">{content}</pre>'
        })

    if sc_items:
        with open("seqcode_dist_mqc.html", "w", encoding="utf-8") as out:
            out.write("<!--\\n")
            out.write("id: seqcode_genome_distribution\\n")
            out.write("section_name: 'SeqCode: Genome Distribution'\\n")
            out.write("description: 'Risultati SeqCode genomeDistribution (scorribili con le frecce).'\\n")
            out.write("-->\\n")
            out.write(build_carousel("seqcode_carousel", "SeqCode Genome Distribution", sc_items))
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
