###############################################################################
# GSE242216 miRNA differential expression: AA vs Healthy Control
# ---------------------------------------------------------------------------
# Dataset      : GEO GSE242216 (NanoString nCounter Human v3 miRNA panel)
# Thresholds   : adj.P < 0.05 AND |log2FC| >= 2 for DE call
#
# Reproducibility
#   R            : >= 4.3
#   Bioconductor : >= 3.18
#   limma        : v3.58+
#   edgeR        : v4.0+
#
# Usage
#   1. Download RCC files from GEO GSE242216 into ./RCC_files/
#   2. Rscript GSE242216_miRNA.R
#   3. Outputs written to ./miRNA_results/
#
###############################################################################

# ---- 0. Dependency check --------------------------------------------------
required_pkgs <- c("limma", "edgeR")
missing_pkgs  <- required_pkgs[!vapply(required_pkgs, requireNamespace,
                                       logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
       "\nInstall with:\n",
       '  if (!requireNamespace("BiocManager", quietly = TRUE)) ',
       'install.packages("BiocManager")\n',
       '  BiocManager::install(c("', paste(missing_pkgs, collapse = '","'),
       '"))')
}

suppressPackageStartupMessages({
  library(limma)
  library(edgeR)
})
set.seed(42)

# ---- USER SETTINGS ---------------------------------------------------------
RCC_DIR <- "RCC_files"
OUT_DIR <- "miRNA_results"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# QC and filter thresholds (reported in Methods)
HK_SF_MIN      <- 0.1
HK_SF_MAX      <- 10
POS_SF_MIN     <- 0.3
POS_SF_MAX     <- 3
EXPR_THRESHOLD <- log2(20)   # mean log2(normalized + 1) filter
LFC_THRESHOLD  <- 2          # |log2FC| threshold for calling DE
ADJP_THRESHOLD <- 0.05

# ---- 1. RCC parser ---------------------------------------------------------
read_rcc <- function(path) {
  con <- if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con))
  lines <- readLines(con, warn = FALSE)
  s <- grep("<Code_Summary>",  lines)
  e <- grep("</Code_Summary>", lines)
  if (length(s) == 0 || length(e) == 0) stop("No Code_Summary block in ", path)
  block <- lines[(s + 1):(e - 1)]
  block <- block[nzchar(block)]
  df <- read.csv(textConnection(block), stringsAsFactors = FALSE)
  df$Count <- suppressWarnings(as.numeric(df$Count))
  df
}

geomean <- function(x) exp(mean(log(pmax(x, 1)), na.rm = TRUE))

# ---- 2. Sample sheet -------------------------------------------------------
rcc_files <- list.files(RCC_DIR, pattern = "\\.RCC(\\.gz)?$", full.names = TRUE)
cat("Found", length(rcc_files), "RCC files\n")
if (length(rcc_files) == 0) stop("No RCCs in ", RCC_DIR)

fn        <- basename(rcc_files)
sample_id <- sub("\\.RCC(\\.gz)?$", "", fn)
parts     <- strsplit(sample_id, "_")
gsm       <- sapply(parts, `[`, 1)
sname     <- sapply(parts, function(x) paste(x[2:(length(x) - 1)], collapse = "_"))
run       <- sapply(parts, function(x) tail(x, 1))
grp       <- ifelse(grepl("^MDS", sname),       "MDS",
             ifelse(grepl("^AA\\.Ctrl", sname), "Control", "AA"))

ss <- data.frame(file = rcc_files, gsm = gsm, sample = sname,
                 run = run, group = grp, stringsAsFactors = FALSE)
cat("\n=== RAW GROUP COUNTS ===\n"); print(table(ss$group))

# Drop MDS — this study is AA vs Control only
ss <- subset(ss, group != "MDS")
ss$group <- factor(ss$group, levels = c("Control", "AA"))
ss$run   <- factor(ss$run)
cat("After MDS drop:", nrow(ss), "samples\n")
print(table(ss$group))

# ---- 3. Read RCCs and build count matrix -----------------------------------
cat("\nReading", nrow(ss), "RCC files...\n")
all_data        <- lapply(ss$file, read_rcc)
names(all_data) <- ss$sample
all_names       <- unique(unlist(lapply(all_data, `[[`, "Name")))

mat <- matrix(NA_real_, nrow = length(all_names), ncol = length(all_data),
              dimnames = list(all_names, ss$sample))
codeclass <- setNames(rep(NA_character_, length(all_names)), all_names)
for (i in seq_along(all_data)) {
  d <- all_data[[i]]
  mat[d$Name, i]    <- d$Count
  codeclass[d$Name] <- d$CodeClass
}
cat("\n=== PROBE CLASS COUNTS ===\n"); print(table(codeclass))

pos_idx  <- which(codeclass == "Positive")
neg_idx  <- which(codeclass == "Negative")
hk_idx   <- which(codeclass == "Housekeeping")
endo_idx <- grep("^Endogenous", codeclass)   # matches Endogenous / Endogenous1 / Endogenous2

cat("Positive:",     length(pos_idx),
    "| Negative:",   length(neg_idx),
    "| Housekeeping:", length(hk_idx),
    "| Endogenous:", length(endo_idx), "\n")

# ---- 4. Standard NanoString normalization ---------------------------------
# (i) Positive-control normalization (geometric mean of ERCC positive probes)
pos_geo <- apply(mat[pos_idx, ], 2, geomean)
pos_sf  <- geomean(pos_geo) / pos_geo
mat_pos <- sweep(mat, 2, pos_sf, "*")

# (ii) Background subtraction (mean of negative-control probes; floor at 1)
neg_mean <- colMeans(mat_pos[neg_idx, ], na.rm = TRUE)
mat_bg   <- sweep(mat_pos, 2, neg_mean, "-")
mat_bg[mat_bg < 1] <- 1

# (iii) Housekeeping normalization (geometric mean of HK probes)
hk_geo   <- apply(mat_bg[hk_idx, ], 2, geomean)
hk_sf    <- geomean(hk_geo) / hk_geo
mat_norm <- sweep(mat_bg, 2, hk_sf, "*")

cat("\n=== SCALE FACTOR SUMMARIES ===\n")
cat("Positive scale factor:\n");    print(summary(pos_sf))
cat("Housekeeping scale factor:\n"); print(summary(hk_sf))

# ---- 5. Sample QC exclusion -----------------------------------------------
qc_flag <- which(hk_sf < HK_SF_MIN | hk_sf > HK_SF_MAX |
                 pos_sf < POS_SF_MIN | pos_sf > POS_SF_MAX)
cat("\n=== QC-FAILED SAMPLES (excluded) ===\n")
cat("N flagged:", length(qc_flag), "\n")
if (length(qc_flag) > 0) {
  qc_tab <- data.frame(sample = ss$sample[qc_flag],
                       group  = ss$group[qc_flag],
                       pos_sf = round(pos_sf[qc_flag], 3),
                       hk_sf  = round(hk_sf[qc_flag], 3))
  print(qc_tab)
  write.csv(qc_tab, file.path(OUT_DIR, "QC_excluded_samples.csv"), row.names = FALSE)
}

# ---- 6. Endogenous miRNA log2 matrix + drop NA probes ---------------------
mirna_mat   <- log2(mat_norm[endo_idx, ] + 1)
n_before_na <- nrow(mirna_mat)
mirna_mat   <- mirna_mat[complete.cases(mirna_mat), ]
cat("\nEndogenous probes:", n_before_na,
    "| after NA-drop:",    nrow(mirna_mat), "\n")

# Apply QC exclusion
if (length(qc_flag) > 0) {
  keep_samp <- setdiff(seq_len(ncol(mirna_mat)), qc_flag)
  mirna_mat <- mirna_mat[, keep_samp]
  ss        <- ss[keep_samp, ]
}
ss$run <- droplevels(ss$run)

cat("\n=== POST-QC COHORT ===\n")
cat("Samples:", ncol(mirna_mat), "\n")
print(table(ss$group))
cat("Run x group:\n"); print(table(ss$run, ss$group))

# ---- 7. Expression filter -------------------------------------------------
keep_mir <- rowMeans(mirna_mat, na.rm = TRUE) >= EXPR_THRESHOLD
cat("\n=== EXPRESSION FILTER (mean log2(norm + 1) >=", EXPR_THRESHOLD, ") ===\n")
cat("miRNAs:", nrow(mirna_mat), "->", sum(keep_mir), "after filter\n")
mirna_mat <- mirna_mat[keep_mir, ]

if (sum(ss$group == "Control") < 2 || sum(ss$group == "AA") < 5)
  stop("Too few samples per group after QC — inspect QC step.")

# ---- 8. Differential expression: limma with ~ run + group ----------------
design <- model.matrix(~ run + group, data = ss)
fit    <- lmFit(mirna_mat, design)
fit    <- eBayes(fit, trend = TRUE)
res    <- topTable(fit, coef = "groupAA", number = Inf, sort.by = "P")
res$miRNA <- rownames(res)
res <- res[, c("miRNA", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]

cat("\n=== AA vs Control | limma ~ run + group ===\n")
cat("Positive logFC = higher in AA\n")
cat("miRNAs tested:", nrow(res), "\n")
cat("adj.P <", ADJP_THRESHOLD, ":",
    sum(res$adj.P.Val < ADJP_THRESHOLD), "\n")
cat("adj.P <", ADJP_THRESHOLD, "AND |logFC| >=", LFC_THRESHOLD, ":",
    sum(res$adj.P.Val < ADJP_THRESHOLD & abs(res$logFC) >= LFC_THRESHOLD),
    "\n")

# Final DE call: FDR AND effect-size threshold
res$DE_call <- ifelse(res$adj.P.Val < ADJP_THRESHOLD &
                        res$logFC  >=  LFC_THRESHOLD, "UP",
               ifelse(res$adj.P.Val < ADJP_THRESHOLD &
                        res$logFC  <= -LFC_THRESHOLD, "DOWN", "NS"))
cat("\nFinal DE call (adjP <", ADJP_THRESHOLD,
    "& |logFC| >=", LFC_THRESHOLD, "):\n")
print(table(res$DE_call))

up_mirs   <- subset(res, DE_call == "UP")
up_mirs   <- up_mirs[order(-up_mirs$logFC), ]
down_mirs <- subset(res, DE_call == "DOWN")
down_mirs <- down_mirs[order(down_mirs$logFC), ]

cat("\n--- UPREGULATED in AA ---\n")
if (nrow(up_mirs)   > 0) print(up_mirs,   row.names = FALSE) else cat("(none)\n")

cat("\n--- DOWNREGULATED in AA ---\n")
if (nrow(down_mirs) > 0) print(down_mirs, row.names = FALSE) else cat("(none)\n")

write.csv(res,       file.path(OUT_DIR, "miRNA_DE_AA_vs_Control_full.csv"), row.names = FALSE)
write.csv(up_mirs,   file.path(OUT_DIR, "miRNA_UP_in_AA.csv"),              row.names = FALSE)
write.csv(down_mirs, file.path(OUT_DIR, "miRNA_DOWN_in_AA.csv"),            row.names = FALSE)

# ---- 9. Session info (for reproducibility) --------------------------------
sink(file.path(OUT_DIR, "sessionInfo.txt"))
print(sessionInfo())
sink()

cat("\n=== DONE ===\n")
cat("Output directory:", normalizePath(OUT_DIR), "\n")
print(list.files(OUT_DIR))
