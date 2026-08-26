###############################################################################
# AA_QC_metrics.R
# -----------------------------------------------------------------------------
# Count-derived QC for the GSE165870 aplastic-anaemia bulk RNA-seq matrix.
#
# Outputs (./qc_results/):
#   qc_per_sample.csv           - library size, gene detection, complexity
#   Fig_PCA_before_filter.png   - PCA on log2(x+1) before low-count filter
#   Fig_PCA_after_filter.png    - PCA on VST after low-count filter
#
# NOTE ON UNAVAILABLE METRICS
# -----------------------------------------------------------------------------
# Mapping / alignment rate and duplication rate CANNOT be derived from a
# GEO-deposited count matrix. Those metrics require the aligner log
# (e.g. STAR Log.final.out) and Picard MarkDuplicates output on the sorted
# BAM, neither of which is contained in GSE165870_Full_length_bulk_counts.txt.
# They are therefore not reported by this script.
#
# What IS reported (from counts alone):
#   - library size per sample (total counts)
#   - number of genes detected (count > 0) per sample
#   - number of genes robustly expressed (count >= 10) per sample
#   - library-complexity proxies: fraction of counts in the top-1 and top-100
#     genes (higher = more skewed transcriptome / possible PCR saturation)
#   - PCA before removing low-count genes  (log2(x+1))
#   - PCA after removing low-count genes    (VST, DESeq2 default)
#
# Author       : Rojina Yasmin
# Affiliation  : Biomedical Genetics Laboratory, Dept. of Zoology,
#                University of Burdwan, West Bengal, India
# License      : MIT
#
# Dependencies:
#   Bioconductor : DESeq2
#   CRAN         : ggplot2, ggrepel
#
#   install.packages("BiocManager")
#   BiocManager::install("DESeq2")
#   install.packages(c("ggplot2", "ggrepel"))
###############################################################################

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
})

set.seed(42)

# -----------------------------------------------------------------------------
# PARAMETERS
# -----------------------------------------------------------------------------
COUNTS_FILE  <- "GSE165870_Full_length_bulk_counts.txt"
OUT_DIR      <- "qc_results"
MIN_ROW_SUM  <- 10       # low-count filter used for the "after-filter" PCA
PCA_NTOP     <- 2000     # top-variable genes for PCA (both before and after)

dir.create(OUT_DIR, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 1. Load count matrix
# -----------------------------------------------------------------------------
raw <- read.table(COUNTS_FILE, header = TRUE, sep = "\t",
                  check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(all(c("GeneName", "GeneSymbol") %in% colnames(raw)))

gene_ids    <- raw$GeneName
sample_cols <- setdiff(colnames(raw), c("GeneName", "GeneSymbol"))
counts_mat  <- as.matrix(raw[, sample_cols])
rownames(counts_mat) <- gene_ids

stopifnot(!any(duplicated(gene_ids)))
stopifnot(!anyNA(counts_mat))
stopifnot(all(counts_mat == round(counts_mat)))
mode(counts_mat) <- "integer"

# -----------------------------------------------------------------------------
# 2. Sample metadata
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

# -----------------------------------------------------------------------------
# 3. Per-sample count-derived QC table
# -----------------------------------------------------------------------------
n_genes_total <- nrow(counts_mat)

qc <- data.frame(
  sample                = colnames(counts_mat),
  condition             = as.character(coldata$condition),
  library_size          = colSums(counts_mat),
  genes_detected_gt0    = colSums(counts_mat > 0),
  frac_genes_detected   = round(colSums(counts_mat > 0) / n_genes_total, 4),
  genes_expressed_ge10  = colSums(counts_mat >= 10),
  median_nonzero_count  = apply(counts_mat, 2,
                                function(x) median(x[x > 0])),
  max_count             = apply(counts_mat, 2, max),
  frac_counts_top1_gene = round(apply(counts_mat, 2,
                                      function(x) max(x) / sum(x)), 4),
  frac_counts_top100    = round(apply(counts_mat, 2, function(x) {
                            sum(sort(x, decreasing = TRUE)[1:100]) / sum(x)
                          }), 4),
  row.names             = NULL
)

cat("\nPer-sample QC (from count matrix only):\n")
print(qc, row.names = FALSE)
write.csv(qc, file.path(OUT_DIR, "qc_per_sample.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 4. PCA BEFORE low-count filtering  (log2(x+1), top-variable genes)
# -----------------------------------------------------------------------------
pre_keep <- rowSums(counts_mat) >= 1     # remove all-zero rows only
pre_mat  <- log2(counts_mat[pre_keep, ] + 1)

pre_vars     <- apply(pre_mat, 1, var)
pre_variable <- pre_mat[pre_vars > 0, ]  # drop zero-variance rows
n_top_pre    <- min(PCA_NTOP, nrow(pre_variable))
pre_top      <- pre_variable[order(-apply(pre_variable, 1, var))[seq_len(n_top_pre)], ]

pca_pre     <- prcomp(t(pre_top), scale. = TRUE)
var_exp_pre <- round(100 * pca_pre$sdev^2 / sum(pca_pre$sdev^2), 1)
pre_df <- data.frame(
  PC1       = pca_pre$x[, 1],
  PC2       = pca_pre$x[, 2],
  sample    = rownames(pca_pre$x),
  condition = coldata[rownames(pca_pre$x), "condition"]
)

p_pre <- ggplot(pre_df, aes(PC1, PC2, color = condition, label = sample)) +
  geom_point(size = 4) +
  geom_text_repel(size = 3.5, show.legend = FALSE, max.overlaps = Inf) +
  scale_color_manual(values = c("Control" = "#377EB8", "AA" = "#E41A1C")) +
  labs(title = sprintf("PCA before filtering  (%d genes, top %d by variance; log2(x+1))",
                       nrow(pre_mat), n_top_pre),
       x = sprintf("PC1 (%s%%)", var_exp_pre[1]),
       y = sprintf("PC2 (%s%%)", var_exp_pre[2])) +
  theme_bw(base_size = 12)
ggsave(file.path(OUT_DIR, "Fig_PCA_before_filter.png"), p_pre,
       width = 6.5, height = 5, dpi = 300)

# -----------------------------------------------------------------------------
# 5. PCA AFTER low-count filtering  (VST, DESeq2 default)
# -----------------------------------------------------------------------------
post_keep   <- rowSums(counts_mat) >= MIN_ROW_SUM
counts_post <- counts_mat[post_keep, ]

dds_qc <- DESeqDataSetFromMatrix(counts_post, coldata, design = ~ condition)
vsd_qc <- vst(dds_qc, blind = TRUE)

post_df <- plotPCA(vsd_qc, intgroup = "condition",
                   ntop = PCA_NTOP, returnData = TRUE)
percentVar_post <- round(100 * attr(post_df, "percentVar"), 1)
if (!"name" %in% names(post_df)) post_df$name <- rownames(post_df)

p_post <- ggplot(post_df, aes(PC1, PC2, color = condition, label = name)) +
  geom_point(size = 4) +
  geom_text_repel(size = 3.5, show.legend = FALSE, max.overlaps = Inf) +
  scale_color_manual(values = c("Control" = "#377EB8", "AA" = "#E41A1C")) +
  labs(title = sprintf("PCA after filtering  (rowSum >= %d: %d genes; VST)",
                       MIN_ROW_SUM, nrow(counts_post)),
       x = sprintf("PC1 (%s%%)", percentVar_post[1]),
       y = sprintf("PC2 (%s%%)", percentVar_post[2])) +
  theme_bw(base_size = 12)
ggsave(file.path(OUT_DIR, "Fig_PCA_after_filter.png"), p_post,
       width = 6.5, height = 5, dpi = 300)

# -----------------------------------------------------------------------------
# Done
# -----------------------------------------------------------------------------
cat("\n=== DONE ===\n")
cat("Wrote to", normalizePath(OUT_DIR), ":\n")
cat("  -", file.path(OUT_DIR, "qc_per_sample.csv"), "\n")
cat("  -", file.path(OUT_DIR, "Fig_PCA_before_filter.png"), "\n")
cat("  -", file.path(OUT_DIR, "Fig_PCA_after_filter.png"), "\n")
