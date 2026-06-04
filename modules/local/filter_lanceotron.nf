process FILTER_LANCEOTRON {
    tag "$meta.id"
    label 'process_low'
    container 'quay.io/biocontainers/python:3.9--1'

    input:
    tuple val(meta), path(peaks)
    val threshold

    output:
    tuple val(meta), path("*.diffbind_ready.bed"), emit: filtered_peaks
    tuple val(meta), path("*.filtered_full.bed") , emit: full_peaks
    tuple val(meta), path("*_filtered_mqc.txt")  , emit: counts_mqc  // <--- NUOVO OUTPUT PER MULTIQC

    script:
    """
    #!/usr/bin/env python3
    import sys

    with open("${peaks}", "r") as f:
        lines = f.readlines()

    out_full = "${meta.id}.filtered_full.bed"
    out_diff = "${meta.id}.diffbind_ready.bed"
    out_mqc  = "${meta.id}_filtered_mqc.txt"

    if not lines:
        open(out_full, "w").close()
        open(out_diff, "w").close()
        # Se non ci sono picchi, stampa 0 per MultiQC
        with open(out_mqc, "w") as mq:
            mq.write("# id: 'lanceotron_filtered'\\n# plot_type: 'generalstats'\\nSample\\tlanceotron_filtered_count\\n${meta.id}\\t0\\n")
        sys.exit(0)

    header = lines[0].strip().split("\\t")
    
    idx = -1
    for i, col in enumerate(header):
        clean_col = col.lower().replace(" ", "_")
        if "overall_peak_score" in clean_col:
            idx = i
            break

    filtered_count = 0  # <--- INIZIALIZZIAMO IL CONTATORE

    with open(out_full, "w") as f_full, open(out_diff, "w") as f_diff:
        if idx == -1:
            f_full.writelines(lines)
            for i, line in enumerate(lines):
                if i == 0 and "chr" not in line.lower():
                    continue
                parts = line.strip().split("\\t")
                if len(parts) >= 3:
                    f_diff.write("\\t".join(parts[:5]) + "\\n")
                    filtered_count += 1
        else:
            f_full.write(lines[0])
            
            for line in lines[1:]:
                parts = line.strip().split("\\t")
                if len(parts) > idx:
                    try:
                        score = float(parts[idx])
                        if score > float(${threshold}):
                            f_full.write(line)
                            filtered_count += 1  # <--- INCREMENTIAMO IL CONTATORE QUANDO SUPERA IL FILTRO
                            
                            chrom = parts[0]
                            start = parts[1]
                            end = parts[2]
                            name = parts[3] if len(parts) > 3 else "${meta.id}_peak"
                            
                            clean_line = chrom + "\\t" + start + "\\t" + end + "\\t" + name + "\\t" + str(score) + "\\n"
                            f_diff.write(clean_line)
                    except ValueError:
                        continue


    with open(out_mqc, "w") as mq:
        mq.write("# id: 'lanceotron_filtered'\\n")
        mq.write("# section_name: 'Lanceotron Filtered Peaks'\\n")
        mq.write("# plot_type: 'generalstats'\\n")
        mq.write("# pconfig:\\n")
        mq.write("#    lanceotron_filtered_count:\\n")
        mq.write("#        title: 'Lanceotron Peaks (Filtered)'\\n")
        mq.write("#        namespace: 'Lanceotron'\\n")
        mq.write("#        desc: 'Numero di picchi che superano la soglia di score impostata'\\n")
        mq.write("#        format: '{:d}'\\n")
        mq.write("Sample\\tlanceotron_filtered_count\\n")
        mq.write(f"${meta.id}\\t{filtered_count}\\n")
    """
}
