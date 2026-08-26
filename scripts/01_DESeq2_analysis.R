###############################################################################
# AA_DESeq2_analysis.R
# -----------------------------------------------------------------------------
# Aplastic Anaemia (AA) bulk RNA-seq differential expression analysis
# Dataset : GSE165870 (Lin-CD34+ HSPCs; 3 healthy controls vs 6 non-SAA cases)
# Method  : DESeq2 (Wald test) with apeglm LFC shrinkage
# Output  : ./deseq2_results/
#             DESeq2_full_results.csv
#             DESeq2_significant_DEGs.csv
#
# -----------------------------------------------------------------------------
# How to run
# -----------------------------------------------------------------------------
#   1. Download the count file from GEO (GSE165870):
#        GSE165870_Full_length_bulk_counts.txt.gz
#      Decompress it and place it in the working directory.
#   2. From an R session in that directory:
#        source("AA_DESeq2_analysis.R")
#      Or run section-by-section in RStudio.
#
# -----------------------------------------------------------------------------
# Dependencies
# -----------------------------------------------------------------------------
#   Bioconductor : DESeq2, apeglm
#
#   install.packages("BiocManager")
#   BiocManager::install(c("DESeq2", "apeglm"))
#
# Tested with R >= 4.3.0.
###############################################################################

suppressPackageStartupMessages({
  library(DESeq2)
})

set.seed(42)  # DESeq2/apeglm are deterministic; kept for downstream reproducibility

# -----------------------------------------------------------------------------
# User-configurable parameters
# -----------------------------------------------------------------------------
COUNTS_FILE  <- "GSE165870_Full_length_bulk_counts.txt"
OUT_DIR      <- "deseq2_results"
MIN_ROW_SUM  <- 10     # pre-filter: keep genes with rowSums >= MIN_ROW_SUM
PADJ_THR     <- 0.05   # BH-adjusted p-value threshold for DEG calling
LFC_THR      <- 2      # |log2 fold change| threshold for DEG calling

dir.create(OUT_DIR, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 1. Load count matrix
# -----------------------------------------------------------------------------
raw <- read.table(COUNTS_FILE, header = TRUE, sep = "\t",
                  check.names = FALSE, stringsAsFactors = FALSE)

cat("Dimensions of input matrix:\n"); print(dim(raw))
cat("Column names:\n"); print(colnames(raw))

# Lock the invariant that the annotation columns are named as expected
stopifnot(all(c("GeneName", "GeneSymbol") %in% colnames(raw)))

# Build count matrix by column NAME (not position) — protects against re-uploads
# with a different column order.
gene_ids    <- raw$GeneName
gene_syms   <- raw$GeneSymbol
sample_cols <- setdiff(colnames(raw), c("GeneName", "GeneSymbol"))
counts_mat  <- as.matrix(raw[, sample_cols])
rownames(counts_mat) <- gene_ids

# Safety checks
stopifnot(length(gene_syms) == nrow(raw))
stopifnot(!any(duplicated(gene_ids)))           # DESeqDataSetFromMatrix requires unique rownames
stopifnot(!anyNA(counts_mat))
stopifnot(all(counts_mat == round(counts_mat))) # confirm true integer counts (not RSEM/Salmon)
mode(counts_mat) <- "integer"

# -----------------------------------------------------------------------------
# 2. Sample metadata (explicit name -> condition mapping)
# -----------------------------------------------------------------------------
sample_labels <- c(
  Ctrl2 = "Control", Ctrl7 = "Control", Ctrl8 = "Control",
  P18   = "AA",      P19   = "AA",      P20   = "AA",
  P21   = "AA",      P22   = "AA",      P23   = "AA"
)
stopifnot(setequal(colnames(counts_mat), names(sample_labels)))

coldata <- data.frame(
  sample    = colnames(counts_mat),
  condition = factor(sample_labels[colnames(counts_mat)],
                     levels = c("Control", "AA")),
  row.names = colnames(counts_mat)
)
cat("\nSample design:\n"); print(coldata)

# Library sizes (Ctrl7 is notably small in raw data — reported in QC)
cat("\nLibrary sizes (total counts per sample):\n")
print(colSums(counts_mat))

# -----------------------------------------------------------------------------
# 3. Pre-filter low-count genes
# -----------------------------------------------------------------------------
keep <- rowSums(counts_mat) >= MIN_ROW_SUM
cat("\nGenes before filtering (rowSums >= ", MIN_ROW_SUM, "): ", nrow(counts_mat),
    "\nGenes after filtering                       : ", sum(keep), "\n", sep = "")
counts_mat <- counts_mat[keep, ]
gene_syms_kept <- gene_syms[keep]
names(gene_syms_kept) <- rownames(counts_mat)

# -----------------------------------------------------------------------------
# 4. DESeq2 differential expression
# -----------------------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = counts_mat,
                              colData   = coldata,
                              design    = ~ condition)

# Print the design matrix so the reader can verify the reference level
cat("\nDesign matrix:\n"); print(model.matrix(~ condition, coldata))

dds <- DESeq(dds)

# Wald test, BH-adjusted p-values, contrast = AA vs Control
res <- results(dds, contrast = c("condition", "AA", "Control"), alpha = PADJ_THR)

# LFC shrinkage (apeglm) for ranking and plotting.
# Note: apeglm shrinks the effect size but reuses the original Wald p-values.
coef_name <- resultsNames(dds)[2]
stopifnot(grepl("AA_vs_Control$", coef_name))
res_shr <- lfcShrink(dds, coef = coef_name, type = "apeglm")

# Build results table with gene symbols
# NOTE: significance is filtered on the SHRUNKEN LFC (|log2FC| >= LFC_THR),
# consistent with using apeglm output for downstream ranking.
res_df <- as.data.frame(res_shr)
res_df$ensembl <- rownames(res_df)
res_df$symbol  <- gene_syms_kept[rownames(res_df)]
res_df <- res_df[, c("ensembl", "symbol", "baseMean", "log2FoldChange",
                     "lfcSE", "pvalue", "padj")]

sig <- subset(res_df,
              !is.na(padj) & padj < PADJ_THR & abs(log2FoldChange) >= LFC_THR)
sig <- sig[order(sig$padj), ]

cat("\nNumber of DEGs (padj<", PADJ_THR,
    ", |log2FC|>=", LFC_THR, "): ", nrow(sig), "\n", sep = "")
cat("  Up   (log2FC >=  ", LFC_THR, "): ", sum(sig$log2FoldChange >=  LFC_THR), "\n", sep = "")
cat("  Down (log2FC <= -", LFC_THR, "): ", sum(sig$log2FoldChange <= -LFC_THR), "\n", sep = "")

# Export full and significant tables
write.csv(res_df, file.path(OUT_DIR, "DESeq2_full_results.csv"),     row.names = FALSE)
write.csv(sig,    file.path(OUT_DIR, "DESeq2_significant_DEGs.csv"), row.names = FALSE)

cat("\n=== DONE ===\n")
cat("Wrote:\n")
cat("  -", file.path(OUT_DIR, "DESeq2_full_results.csv"), "\n")
cat("  -", file.path(OUT_DIR, "DESeq2_significant_DEGs.csv"), "\n")
