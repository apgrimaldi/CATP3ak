# CATP3ak

### ChiP-Seq and ATAC-Seq callpeakers

A Nextflow DSL2 pipeline for comprehensive ATAC-seq and ChIP-seq data analysis.

CATP3ak was developed to improve peak identification by integrating multiple complementary peak-calling approaches within a single reproducible workflow. Starting from raw sequencing data, the pipeline performs quality control, preprocessing, alignment, duplicate removal, blacklist filtering, peak calling, annotation, differential binding analysis, signal profiling, and reporting.

By combining statistical, deep-learning, and Hidden Markov Model-based approaches, CATP3ak enables a comprehensive evaluation of chromatin accessibility and protein–DNA interaction data while facilitating downstream biological interpretation.


[![Nextflow](https://img.shields.io/badge/Nextflow-DSL2-brightgreen)](https://www.nextflow.io/)
[![Docker](https://img.shields.io/badge/Container-Docker-blue)](https://www.docker.com/)

---

## Overview

CATP3ak is a modular workflow designed for both ChIP-seq and ATAC-seq experiments.
The pipeline automatically:

Detects **Single-End** and **Paired-End** libraries
Identifies control samples for ChIP-seq analyses
Supports custom reference genomes
Generates browser-ready signal tracks
Executes multiple peak-calling strategies
Produces publication-ready quality-control reports

---
Features
ChIP-seq and ATAC-seq support
Automatic Single-End / Paired-End detection
Automatic control sample identification
Quality control with FastQC
Adapter trimming with Trim Galore
Alignment using Bowtie2
BAM processing with SAMtools
Duplicate removal with Picard
Blacklist filtering
BigWig signal track generation
Multi-strategy peak calling:
MACS3 (statistical modeling)
Lanceotron (deep learning)
OmniPeak (Hidden Markov Models)
FRiP score calculation
Peak annotation with HOMER
Differential binding analysis with DiffBind
Signal profiling with Profileplyr
MultiQC reporting
Docker support
AWS/S3 compatibility through nf-amazon

---

## Workflow

```text
FASTQ
 ↓
FastQC
 ↓
Trim Galore
 ↓
Bowtie2 Alignment
 ↓
SAMtools Processing
 ↓
Picard MarkDuplicates
 ↓
Blacklist Filtering
 ↓
deepTools QC + BigWig Generation
 ↓
 ├── MACS3 Peak Calling
 └── Lanceotron Peak Calling
          ↓
       FRiP Score
          ↓
     HOMER Annotation
          ↓
        DiffBind
          ↓
      Profileplyr
          ↓
         MultiQC
```

---

## Requirements

* Nextflow ≥ 25.x
* Docker (recommended) or Singularity/Apptainer
* Linux environment

Verify installation:

```bash
nextflow -version
docker --version
```

---

## Quick Start

Run the pipeline directly from GitHub:

```bash
nextflow run apgrimaldi/CATP3ak \
    -latest \
    -profile docker \
    --input samplesheet.csv \
    --protocol chip \
    --genome GRCh38 \
    --outdir results
```

---

## Main Parameters

| Parameter                | Description                                       |
| ------------------------ | ------------------------------------------------- |
| `--input`                | Input samplesheet                                 |
| `--protocol`             | Analysis type (`chip` or `atac`)                  |
| `--genome`               | Genome identifier                                 |
| `--outdir`               | Output directory                                  |
| `--fragment_size`        | Fragment size used for single-end analyses        |
| `--single_end`           | Force single-end processing                       |
| `--lanceotron_threshold` | Minimum Lanceotron score retained after filtering |
| `--skip_homer`           | Skip HOMER annotation                             |
| `--skip_diffbind`        | Skip DiffBind analysis                            |
| `--skip_profileplyr`     | Skip Profileplyr analysis                         |

---

## Custom Genome Support

CATP3ak supports custom reference genomes.

### Parameters

| Parameter         | Description                           |
| ----------------- | ------------------------------------- |
| `--fasta_file`    | Reference genome FASTA file           |
| `--gtf_file`      | Gene annotation GTF file              |
| `--macs_gsize`    | Effective genome size for MACS3       |
| `--blacklist`     | BED file containing blacklist regions |
| `--bowtie2_index` | Pre-built Bowtie2 index               |

### Example

```bash
nextflow run apgrimaldi/CATP3ak \
    -profile docker \
    --protocol chip \
    --input samplesheet.csv \
    --fasta_file reference.fasta \
    --gtf_file annotation.gtf \
    --macs_gsize 2.7e9 \
    --blacklist blacklist.bed \
    --outdir results \
    -resume
```

---

## Input Samplesheet

### ChIP-seq Example

```csv
sample,fastq_1,fastq_2,antibody,control,is_control
IP_H3K27ac_1,data/IP_H3K27ac_1.fastq.gz,,H3K27ac,Input_1,false
Input_1,data/Input_1.fastq.gz,,, ,true
```

### ATAC-seq Example

```csv
sample,fastq_1,fastq_2
ATAC_1,data/ATAC_1_R1.fastq.gz,data/ATAC_1_R2.fastq.gz
ATAC_2,data/ATAC_2_R1.fastq.gz,data/ATAC_2_R2.fastq.gz
```

### Samplesheet Columns

| Column       | Description                                         |
| ------------ | --------------------------------------------------- |
| `sample`     | Unique sample identifier                            |
| `fastq_1`    | Read 1 FASTQ file                                   |
| `fastq_2`    | Read 2 FASTQ file (leave empty for Single-End data) |
| `antibody`   | Antibody name (ChIP-seq only)                       |
| `control`    | Matching control sample                             |
| `is_control` | Optional control flag (`true` or `false`)           |

---

## Automatic Control Detection

For ChIP-seq analyses, CATP3ak automatically identifies control samples using one or more of the following criteria:

1. The sample is referenced in the `control` column.
2. The antibody is specified as `IgG`.
3. The `is_control` column is set to `true`.

No control detection is performed in ATAC-seq mode.

---

## Output Structure

```text
results/
├── 00_genome_index/
├── 01_fastqc/
├── 02_trimmed/
├── 03_aligned/
│   ├── raw_bam/
│   ├── sorted_bam/
│   ├── indexed_sorted_bam/
│   └── stats/
├── 04_duplicates_removed/
├── 05_final_filtered_bam/
├── 06_bigwig/
│   └── qc_fingerprint/
├── 07_lanceotron/
│   ├── unfiltered/
│   ├── filtered/
│   └── bigwig_res1/
├── 08_peaks_macs3/
│   ├── narrow/
│   ├── broad/
│   └── frip_stats/
├── 09_annotation/
│   ├── macs/
│   └── lanceotron/
├── 10_diffbind/
│   ├── macs/
│   └── lanceotron/
├── 11_profileplyr/
│   ├── macs/
│   └── lanceotron/
└── 12_MultiQC_Report/
```

---

## Generated Outputs

CATP3ak produces:

* Quality-control reports
* Filtered BAM files
* BigWig tracks for genome browsers
* MACS3 peaks (narrow and broad)
* Lanceotron peaks (raw and filtered)
* FRiP statistics
* HOMER annotations
* DiffBind differential binding results
* Profileplyr signal profiling reports
* MultiQC summary report

---

## Reproducibility

All software dependencies are executed within containers, ensuring reproducible analyses across different computational environments.

Supported execution environments:

* Docker
* Singularity / Apptainer
* AWS-compatible infrastructures

---

## Author

**Annapaola Grimaldi**

Biological Sciences

GitHub: https://github.com/apgrimaldi

---

CATP3ak was developed to provide a reproducible and user-friendly framework for chromatin accessibility and protein–DNA interaction studies.
