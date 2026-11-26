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
    # --- 1) DeepTools TSS Heatmaps (*.png) con carosello a frecce ---
    shopt -s nullglob
    dt_pngs=(heatmaps_in/*.png)
    if [ \${#dt_pngs[@]} -gt 0 ]; then
      cat << 'EOF' > deeptools_heatmaps_mqc.html
<!--
id: deeptools_tss_heatmaps
section_name: 'DeepTools: TSS Heatmaps'
description: 'Heatmap del segnale attorno ai TSS (scorribili con le frecce).'
-->
<div style="border:1px solid #ddd; border-radius:8px; padding:16px; margin:12px 0;">
  <div style="display:flex; justify-content:center; align-items:center; gap:16px; margin-bottom:12px;">
    <button type="button" onclick="move_dt_carousel(-1)" style="padding:6px 14px; cursor:pointer; font-size:16px;">&#9664; Prev</button>
    <span id="dt_carousel_counter" style="font-weight:600;">1 / \${#dt_pngs[@]}</span>
    <button type="button" onclick="move_dt_carousel(1)" style="padding:6px 14px; cursor:pointer; font-size:16px;">Next &#9654;</button>
  </div>
EOF

      idx=0
      total=\${#dt_pngs[@]}
      for p in "\${dt_pngs[@]}"; do
        fname=\$(basename "\$p")
        b64=\$(base64 -w 0 "\$p")
        disp="none"
        [ \$idx -eq 0 ] && disp="block"
        num=\$((idx + 1))
        cat << EOF >> deeptools_heatmaps_mqc.html
  <div class="dt_carousel_slide" style="display:\${disp}; text-align:center;">
    <div style="font-weight:600; margin-bottom:8px;">\${fname} (\${num}/\${total})</div>
    <img src="data:image/png;base64,\${b64}" style="max-width:100%; max-height:720px; height:auto;" />
  </div>
EOF
        idx=\$((idx + 1))
      done

      cat << 'EOF' >> deeptools_heatmaps_mqc.html
</div>
<script>
(function() {
  let idx = 0;
  const slides = document.getElementsByClassName("dt_carousel_slide");
  const counter = document.getElementById("dt_carousel_counter");
  window.move_dt_carousel = function(step) {
    if (!slides.length) return;
    slides[idx].style.display = "none";
    idx = (idx + step + slides.length) % slides.length;
    slides[idx].style.display = "block";
    if (counter) counter.textContent = (idx + 1) + " / " + slides.length;
  };
})();
</script>
EOF
    fi

    # --- 2) SeqCode Genome Distribution (*.pdf) con carosello a frecce ---
    sc_pdfs=(seqcode_in/*.pdf)
    if [ \${#sc_pdfs[@]} -gt 0 ]; then
      cat << 'EOF' > seqcode_dist_mqc.html
<!--
id: seqcode_genome_distribution
section_name: 'SeqCode: Genome Distribution'
description: 'Report PDF di SeqCode genomeDistribution (scorribili con le frecce).'
-->
<div style="border:1px solid #ddd; border-radius:8px; padding:16px; margin:12px 0;">
  <div style="display:flex; justify-content:center; align-items:center; gap:16px; margin-bottom:12px;">
    <button type="button" onclick="move_sc_carousel(-1)" style="padding:6px 14px; cursor:pointer; font-size:16px;">&#9664; Prev</button>
    <span id="sc_carousel_counter" style="font-weight:600;">1 / \${#sc_pdfs[@]}</span>
    <button type="button" onclick="move_sc_carousel(1)" style="padding:6px 14px; cursor:pointer; font-size:16px;">Next &#9654;</button>
  </div>
EOF

      idx=0
      total=\${#sc_pdfs[@]}
      for pdf in "\${sc_pdfs[@]}"; do
        fname=\$(basename "\$pdf")
        b64=\$(base64 -w 0 "\$pdf")
        disp="none"
        [ \$idx -eq 0 ] && disp="block"
        num=\$((idx + 1))
        cat << EOF >> seqcode_dist_mqc.html
  <div class="sc_carousel_slide" style="display:\${disp}; text-align:center;">
    <div style="font-weight:600; margin-bottom:8px;">\${fname} (\${num}/\${total})</div>
    <embed src="data:application/pdf;base64,\${b64}" type="application/pdf" style="width:100%; height:680px; border:1px solid #ccc;" />
    <div style="margin-top:8px; font-size:13px;">
      <a download="\${fname}" href="data:application/pdf;base64,\${b64}">Scarica \${fname}</a>
    </div>
  </div>
EOF
        idx=\$((idx + 1))
      done

      cat << 'EOF' >> seqcode_dist_mqc.html
</div>
<script>
(function() {
  let idx = 0;
  const slides = document.getElementsByClassName("sc_carousel_slide");
  const counter = document.getElementById("sc_carousel_counter");
  window.move_sc_carousel = function(step) {
    if (!slides.length) return;
    slides[idx].style.display = "none";
    idx = (idx + step + slides.length) % slides.length;
    slides[idx].style.display = "block";
    if (counter) counter.textContent = (idx + 1) + " / " + slides.length;
  };
})();
</script>
EOF
    fi

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
