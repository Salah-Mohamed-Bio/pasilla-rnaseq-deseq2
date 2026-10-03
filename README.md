# RNA-Seq Reference-Based Analysis — *Pasilla* Gene Knockdown in *Drosophila melanogaster*

A complete, reproducible RNA-Seq differential expression pipeline, implemented entirely via the command line and R, based on the [Galaxy Training Network tutorial](https://training.galaxyproject.org/training-material/topics/transcriptomics/tutorials/ref-based/tutorial.html).

## Overview

This project re-analyzes a subset of the data from [Brooks et al. 2011](https://doi.org/10.1101/gr.108662.110), who studied the effect of depleting the *Pasilla* gene (the *Drosophila* homolog of the mammalian splicing regulators Nova-1/Nova-2) using RNA interference (RNAi).

**Biological question:** Which genes and pathways change expression when *Pasilla* is knocked down?

Two conditions are compared:
- **Untreated** — normal cells
- **Treated** — *Pasilla* gene knocked down via RNAi

## Pipeline

| Step | Tool | Purpose |
|---|---|---|
| 1. Data download | `wget` (Zenodo) | Raw paired-end FASTQ files |
| 2. Quality control | Falco + MultiQC | Per-base quality, GC content, duplication |
| 3. Trimming | Cutadapt | Remove low-quality bases / short reads |
| 4. Mapping | STAR | Splice-aware alignment to the *D. melanogaster* genome (dm6) |
| 5. Read counting | featureCounts | Reads per gene, per sample |
| 6. Differential expression | DESeq2 (R) | Identify significantly up/down-regulated genes |
| 7. Annotation | AnnotationDbi + org.Dm.eg.db | Map FlyBase IDs to gene symbols/names |
| 8. Visualization | plotMA, EnhancedVolcano, pheatmap | MA plot, Volcano plot, Heatmap |
| 9. Functional enrichment | clusterProfiler | GO Biological Process enrichment |

## Repository Structure

```
pasilla-rnaseq-deseq2/
├── README.md
├── deseq2_analysis.R            # DESeq2 + annotation + plots + GO enrichment
└── results/
    ├── deseq2_full_results.csv
    ├── significant_genes.csv
    ├── go_enrichment_results.csv
    ├── MA_plot.png
    ├── volcano_plot.png
    ├── heatmap.png
    └── go_dotplot.png
```

> **Note:** Only the R analysis script is included in this repository. The upstream command-line steps (QC, trimming, mapping, read counting) were run locally following the commands documented below, but are not committed as standalone scripts here. Raw FASTQ files, BAM files, and the reference genome/GTF are **not** included due to their size — see [Data Sources](#data-sources) for download links.

## Requirements

### Command-line tools (via conda/mamba)

```bash
mamba create -n rnaseq -c bioconda -c conda-forge \
    falco multiqc cutadapt star subread samtools -y
mamba activate rnaseq
```

### R packages

```r
install.packages("BiocManager")
BiocManager::install(c(
    "DESeq2", "AnnotationDbi", "org.Dm.eg.db",
    "clusterProfiler", "EnhancedVolcano"
))
install.packages("pheatmap")
```

> System libraries `libcairo2-dev` and `libuv1-dev` are required on Ubuntu/Debian for some R dependencies (`gdtools`, `fs`):
> ```bash
> sudo apt install -y libcairo2-dev libuv1-dev
> ```

## Data Sources

| File | Source |
|---|---|
| FASTQ reads (subsampled) | [Zenodo record 6457007](https://zenodo.org/record/6457007) |
| FASTQ reads (full, 7 samples) | Same Zenodo record |
| Reference genome (dm6) | [UCSC Genome Browser](https://genome.ucsc.edu/) |
| Gene annotation (GTF) | [Zenodo record 6457007](https://zenodo.org/record/6457007/files/Drosophila_melanogaster.BDGP6.32.109_UCSC.gtf.gz) |
| Pre-computed featureCounts tables (7 samples) | Same Zenodo record |

## How to Reproduce

### 1. Clone this repo and set up the environment

```bash
git clone https://github.com/<your-username>/pasilla-rnaseq-deseq2.git
cd pasilla-rnaseq-deseq2

mamba create -n rnaseq -c bioconda -c conda-forge \
    falco multiqc cutadapt star subread samtools -y
mamba activate rnaseq
```

### 2. Download the data

```bash
wget https://zenodo.org/record/6457007/files/GSM461177_1_subsampled.fastqsanger
wget https://zenodo.org/record/6457007/files/GSM461177_2_subsampled.fastqsanger
wget https://zenodo.org/record/6457007/files/GSM461180_1_subsampled.fastqsanger
wget https://zenodo.org/record/6457007/files/GSM461180_2_subsampled.fastqsanger
```

### 3. Quality control

```bash
mkdir -p qc
falco GSM461177_1_subsampled.fastqsanger -o qc/
falco GSM461177_2_subsampled.fastqsanger -o qc/
falco GSM461180_1_subsampled.fastqsanger -o qc/
falco GSM461180_2_subsampled.fastqsanger -o qc/
multiqc qc/ -o qc/multiqc_report
```

### 4. Trimming

```bash
cutadapt -q 20 -m 20 \
    -o GSM461177_1_trimmed.fastq -p GSM461177_2_trimmed.fastq \
    GSM461177_1_subsampled.fastqsanger GSM461177_2_subsampled.fastqsanger

cutadapt -q 20 -m 20 \
    -o GSM461180_1_trimmed.fastq -p GSM461180_2_trimmed.fastq \
    GSM461180_1_subsampled.fastqsanger GSM461180_2_subsampled.fastqsanger
```

### 5. Mapping (STAR)

```bash
STAR --runMode genomeGenerate --genomeDir star_index \
    --genomeFastaFiles dm6.fasta \
    --sjdbGTFfile Drosophila_melanogaster.BDGP6.32.109_UCSC.gtf \
    --sjdbOverhang 36

STAR --genomeDir star_index \
    --readFilesIn GSM461177_1_trimmed.fastq GSM461177_2_trimmed.fastq \
    --outSAMtype BAM SortedByCoordinate \
    --quantMode GeneCounts \
    --outFileNamePrefix mapping/untreated_
```

### 6. Read counting (featureCounts)

```bash
featureCounts -a Drosophila_melanogaster.gtf \
    -o count/gene_counts.txt -p -T 4 \
    mapping/treated.sorted.bam mapping/untreated.sorted.bam
```

### 7. Differential expression + enrichment analysis (R)

```bash
Rscript deseq2_analysis.R
```

## Key Results

- **18,589** genes with non-zero counts were tested.
- **966** genes were significantly differentially expressed (adjusted p-value < 0.1):
  - **507** up-regulated
  - **459** down-regulated
- **208** genes passed a stricter significance threshold (padj < 0.05) and were used for functional enrichment.
- GO Biological Process enrichment identified **7 pathways**, most notably:
  - *defense response to Gram-positive bacterium* (strongest signal, FoldEnrichment ≈ 7.7)
  - *dopamine metabolic process*
  - *catecholamine metabolic process*

> **Caveat:** Enrichment results are based on raw p-values (not adjusted) due to the relatively small number of significant genes (208) in this subset analysis. They should be treated as exploratory rather than conclusive. A full re-analysis using all 7 samples from the original study would yield more statistically robust results.

## Acknowledgments

Pipeline and tutorial structure adapted from the [Galaxy Training Network](https://training.galaxyproject.org/) (Batut et al., Freeberg et al., and contributors). Original study: Brooks AN, Yang L, Duff MO, et al. (2011) *Conservation of an RNA regulatory map between Drosophila and mammals*. Genome Research.

## License

This repository's code is shared under the MIT License. The underlying data is subject to the terms of the original Zenodo deposition (CC-BY).
