###############################################################################
# 02_QC_metrics.R
# -----------------------------------------------------------------------------
# Count-derived QC + before/after low-count-filter PCA for the GSE165870
# aplastic-anaemia bulk RNA-seq matrix.
#
# HOW TO RUN LOCALLY
# -----------------------------------------------------------------------------
# 1. Place this script in the folder containing
#      GSE165870_Full_length_bulk_counts.txt
# 2. setwd() to that folder, or use RStudio: Session > Set Working Directory >
#    To Source File Location.
# 3. source("02_QC_metrics.R")
#    First run installs any missing packages automatically.
# 4. Results are written to ./qc_results/ (created if it does not exist).
#
# Outputs (./qc_results/):
#   qc_per_sample.csv           - library size, gene detection, complexity
#   qc_library_size_summary.txt - min / median / max lib size (for Methods)
#   Fig_PCA_before_filter.png   - PCA on rowSums >= 1 matrix (VST)
#   Fig_PCA_after_filter.png    - PCA on rowSums >= MIN_ROW_SUM matrix (VST)
#
# Dependencies
# -----------------------------------------------------------------------------
#   Bioconductor : DESeq2
#   CRAN         : ggplot2, ggrepel
#
###############################################################################

# -----------------------------------------------------------------------------
# 0. Bootstrap: install any missing packages, then load
# -----------------------------------------------------------------------------
.install_if_missing <- function(cran = character(), bioc = character()) {
  need_cran <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
  if (length(need_cran)) {
    message("Installing CRAN packages: ", paste(need_cran, collapse = ", "))
    install.packages(need_cran, repos = "https://cloud.r-project.org")
  }
  need_bioc <- bioc[!vapply(bioc, requireNamespace, logical(1), quietly = TRUE)]
  if (length(need_bioc)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    }
    message("Installing Bioconductor packages: ",
            paste(need_bioc, collapse = ", "))
    BiocManager::install(need_bioc, update = FALSE, ask = FALSE)
  }
}
.install_if_missing(cran = c("ggplot2", "ggrepel"),
                    bioc = c("DESeq2"))

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
MIN_ROW_SUM  <- 10       # low-count filter used in Methods
PCA_NTOP     <- 2000     # top-variable genes for both PCA panels

# -----------------------------------------------------------------------------
# Pre-flight check
# -----------------------------------------------------------------------------
if (!file.exists(COUNTS_FILE)) {
  stop("\nCount matrix not found:\n  ",
       normalizePath(COUNTS_FILE, mustWork = FALSE),
       "\nWorking directory is: ", getwd(),
       "\nMove the file here, or setwd() to its folder.", call. = FALSE)
}
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

n_before_filter <- nrow(counts_mat)
cat("Loaded counts:", n_before_filter, "genes x",
    ncol(counts_mat), "samples\n")

# -----------------------------------------------------------------------------
# 2. Sample metadata
# -----------------------------------------------------------------------------
sample_labels <- c(
  Ctrl2 = "Control", Ctrl7 = "Control", Ctrl8 = "Control",
  P18   = "AA",      P19   = "AA",      P20   = "AA",
  P21   = "AA",      P22   = "AA",      P23   = "AA"
)
if (!setequal(colnames(counts_mat), names(sample_labels))) {
  stop("\nSample names in the count matrix do not match the expected set.\n",
       "Expected: ", paste(names(sample_labels), collapse = ", "), "\n",
       "Found   : ", paste(colnames(counts_mat), collapse = ", "),
       call. = FALSE)
}

coldata <- data.frame(
  sample    = colnames(counts_mat),
  condition = factor(sample_labels[colnames(counts_mat)],
                     levels = c("Control", "AA")),
  row.names = colnames(counts_mat)
)

lib_size <- colSums(counts_mat)

# -----------------------------------------------------------------------------
# 3. Per-sample count-derived QC table  (computed on the RAW matrix)
# -----------------------------------------------------------------------------
qc <- data.frame(
  sample                = colnames(counts_mat),
  condition             = as.character(coldata$condition),
  library_size          = lib_size,
  genes_detected_gt0    = colSums(counts_mat > 0),
  frac_genes_detected   = round(colSums(counts_mat > 0) / n_before_filter, 4),
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

# Library-size summary — the numbers quoted in Methods
libsize_lines <- c(
  "Library size summary (raw counts, per sample) - GSE165870",
  strrep("-", 60),
  sprintf("N samples : %d", length(lib_size)),
  sprintf("Minimum   : %d  (%s)",   min(lib_size),    names(which.min(lib_size))),
  sprintf("Median    : %d",         as.integer(median(lib_size))),
  sprintf("Maximum   : %d  (%s)",   max(lib_size),    names(which.max(lib_size))),
  sprintf("Range     : %.1f-fold",  max(lib_size) / min(lib_size))
)
writeLines(libsize_lines, file.path(OUT_DIR, "qc_library_size_summary.txt"))
cat("\n"); cat(libsize_lines, sep = "\n"); cat("\n")

# -----------------------------------------------------------------------------
# Helper: run the DE pipeline (DESeq -> vst blind=FALSE) and build a PCA plot
# -----------------------------------------------------------------------------
# blind = FALSE means the variance-stabilizing transformation uses the
# same size factors and dispersion trend that DESeq() estimated for the DE
# model. This aligns the PCA with the Methods-described pipeline.
vst_pca_plot <- function(counts_sub, coldata_sub, ntop, title) {
  dds <- DESeqDataSetFromMatrix(counts_sub, coldata_sub, design = ~ condition)
  dds <- DESeq(dds)                       # size factors + dispersions + Wald
  vsd <- vst(dds, blind = FALSE)          # design-informed VST
  d   <- plotPCA(vsd, intgroup = "condition", ntop = ntop, returnData = TRUE)
  pv  <- round(100 * attr(d, "percentVar"), 1)
  if (!"name" %in% names(d)) d$name <- rownames(d)

  ggplot(d, aes(PC1, PC2, color = condition, label = name)) +
    geom_point(size = 4) +
    geom_text_repel(size = 3.5, show.legend = FALSE, max.overlaps = Inf) +
    scale_color_manual(values = c("Control" = "#377EB8", "AA" = "#E41A1C")) +
    labs(title = title,
         x = sprintf("PC1 (%s%%)", pv[1]),
         y = sprintf("PC2 (%s%%)", pv[2])) +
    theme_bw(base_size = 12)
}

# -----------------------------------------------------------------------------
# 4. PCA BEFORE low-count filtering  (rowSums >= 1, only all-zero rows removed)
# -----------------------------------------------------------------------------
pre_keep   <- rowSums(counts_mat) >= 1
counts_pre <- counts_mat[pre_keep, ]
cat("\nBefore filter (rowSum >= 1): ", nrow(counts_pre), " of ",
    n_before_filter, " genes retained.\n", sep = "")

cat("Running DESeq() + VST for 'before filter' PCA...\n")
p_pre <- vst_pca_plot(counts_pre, coldata, ntop = PCA_NTOP,
                      title = sprintf("before filtering (rowSum \u2265 1; %d genes)",
                                      nrow(counts_pre)))
ggsave(file.path(OUT_DIR, "Fig_PCA_before_filter.png"), p_pre,
       width = 6.5, height = 5, dpi = 300)

# -----------------------------------------------------------------------------
# 5. PCA AFTER low-count filtering  (rowSums >= MIN_ROW_SUM, matches Methods)
# -----------------------------------------------------------------------------
post_keep   <- rowSums(counts_mat) >= MIN_ROW_SUM
counts_post <- counts_mat[post_keep, ]
cat("\nAfter filter (rowSum >= ", MIN_ROW_SUM, "): ", nrow(counts_post),
    " of ", n_before_filter, " genes retained.\n", sep = "")

cat("Running DESeq() + VST for 'after filter' PCA...\n")
p_post <- vst_pca_plot(counts_post, coldata, ntop = PCA_NTOP,
                       title = sprintf("after filtering (rowSum \u2265 %d; %d genes)",
                                       MIN_ROW_SUM, nrow(counts_post)))
ggsave(file.path(OUT_DIR, "Fig_PCA_after_filter.png"), p_post,
       width = 6.5, height = 5, dpi = 300)

# -----------------------------------------------------------------------------
# Done
# -----------------------------------------------------------------------------
cat("\n=== DONE ===\n")
cat("Wrote to", normalizePath(OUT_DIR), ":\n")
for (f in c("qc_per_sample.csv",
            "qc_library_size_summary.txt",
            "Fig_PCA_before_filter.png",
            "Fig_PCA_after_filter.png")) {
  cat("  -", file.path(OUT_DIR, f), "\n")
}
