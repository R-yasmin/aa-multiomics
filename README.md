# aa-multiomics

## Overview

This repository contains the analysis pipeline for an *in silico* secondary analysis of two publicly available GEO datasets, integrating bulk RNA-seq of hematopoietic stem/progenitor cells with NanoString miRNA profiling to characterize transcriptomic and microRNA signatures in aplastic anaemia (AA).

## Datasets

| Accession | Platform | Cohort |
|---|---|---|
| [GSE165870](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165870) | Bulk RNA-seq, Lin⁻CD34⁺ HSPCs | 6 AA vs 3 healthy controls |
| [GSE242216](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE242216) | NanoString nCounter Human v3 miRNA panel | AA vs healthy controls (MDS excluded) |

Both datasets are publicly available on GEO.

## Analysis pipeline

1. **Differential gene expression** (`01_DESeq2_analysis.R`) — DESeq2 (Wald test) with apeglm log fold-change shrinkage on GSE165870 counts. DEGs called at adj. p < 0.05 and |log2FC| ≥ 2.
2. **QC of bulk RNA seq** (`02_QC_metrics.R`) — Per-sample library size, and PCA before and after low-count filtering for GSE165870.
3. **miRNA differential expression** (`03_NanoString_miRNA.R`) — Standard NanoString normalization (positive-control, background subtraction, housekeeping) followed by limma with cartridge run as a fixed-effect covariate. DE call at adj. p < 0.05 and |log2FC| ≥ 2.
4. **miRNA–mRNA integration** (`04_multimir_integration.R`) — multiMiR queries across 8 prediction (TargetScan, miRDB, miRanda, DIANA-microT, ElMMo, PITA, PicTar, MicroCosm) and 4 validation (miRTarBase, TarBase, miRecords, miR2Disease) databases. High-confidence pairs = experimentally validated OR predicted by ≥ 2 databases.

## Repository structure
aa-multiomics/
├── scripts/
│ ├── 01_DESeq2_analysis.R
│ ├── 02_QC_metrics.R
│ ├── 03_NanoString_miRNA.R
│ └── 04_multimir_integration.R
└── results/
├── deseq2_results/ # Full and significant DEG tables
├── qc_results/ # QC metrics and PCA figures
├── miRNA_results/ # miRNA DE tables
└── integration_results/ # miRNA–mRNA pair tables and summaries

## Reproducing the analysis

1. Download `GSE165870_Full_length_bulk_counts.txt` from GEO and place it in `scripts/`.
2. Download RCC files from GSE242216 and place them in `scripts/RCC_files/`.
3. In R, from the repository root:

```r
setwd("scripts")
source("01_DESeq2_analysis.R")
source("02_QC_metrics.R")
source("03_NanoString_miRNA.R")
source("04_multimir_integration.R")
```

## Requirements

- R ≥ 4.3.3
- Bioconductor 3.18
- CRAN: `ggplot2`, `ggrepel`, `dplyr`, `tidyr`
- Bioconductor: `DESeq2`, `apeglm`, `limma`, `edgeR`, `multiMiR`

Installation:

```r
install.packages("BiocManager")
BiocManager::install(c("DESeq2", "apeglm", "limma", "edgeR", "multiMiR"))
install.packages(c("ggplot2", "ggrepel", "dplyr", "tidyr"))
```

## Contact

Author Rojina Yasmin (rojinayasmin1996@gmail.com).
