###############################################################################
# miRNA-mRNA integration analysis
#
# Intersects:
#   (a) downregulated miRNAs from GSE242216 limma analysis
#   (b) DEGs from GSE165870 DESeq2 analysis
#   (c) miRNA-target relationships from multiMiR restricted to the 11
#       databases named in the manuscript methodology:
#         Predicted: TargetScan, miRDB, miRanda, DIANA-microT, ElMMo, PITA,
#                    PicTar, MicroCosm
#         Validated: miRTarBase, TarBase, miRecords, miR2Disease
#
# High confidence criterion:
#   experimentally validated (multiMiR type == "validated")
#   OR predicted by >= 2 independent prediction databases
#
# Output (in integration_results/):
#   miRNA-target_all_pairs.csv         all pairs with DB evidence
#   miRNA-target_predicted.csv         predicted pairs only
#   miRNA-target_validated.csv         experimentally validated pairs
#   miRNA-target_high_confidence.csv   validated OR >=2 predicted DBs
#   summary_per_miRNA.csv              per-miRNA target counts
#   summary_per_target.csv             per-target miRNA counts
#
# Setup:
#   renv::install("bioc::multiMiR")
###############################################################################

suppressPackageStartupMessages({
  library(multiMiR)
  library(dplyr)
  library(tidyr)
})
set.seed(42)

cat("multiMiR version:", as.character(packageVersion("multiMiR")), "\n")

# ---- USER SETTINGS ---------------------------------------------------------
DEG_FILE <- "DESeq2_significant_DEGs.csv"        # mRNA DEGs from DESeq2
DEM_FILE <- "miRNA_results/miRNA_DOWN_in_AA.csv" # downregulated miRNAs

OUT_DIR  <- "integration_results"
dir.create(OUT_DIR, showWarnings = FALSE)

# Thresholds
LFC_UP       <- 2            # log2FC >= 2, matches DESeq2 manuscript cutoff
PADJ_THR     <- 0.05
MIN_DBS_PRED <- 2            # min predicted databases for "high confidence"

# The 11 databases named in the manuscript methodology
PRED_DBS <- c("targetscan","mirdb","miranda","diana_microt",
              "elmmo","pita","pictar","microcosm")
VAL_DBS  <- c("mirtarbase","tarbase","mirecords","mir2disease")
METHOD_DBS <- c(PRED_DBS, VAL_DBS)

# ---- 1. Load DEG and DEM lists ---------------------------------------------
cat("Loading DEGs from", DEG_FILE, "...\n")
deg <- read.csv(DEG_FILE, stringsAsFactors = FALSE)
cat("Total DEGs in file:", nrow(deg), "\n")

deg_up <- subset(deg, !is.na(padj) & padj < PADJ_THR &
                      log2FoldChange >= LFC_UP &
                      !is.na(symbol) & symbol != "")
deg_up <- deg_up[order(-deg_up$log2FoldChange), ]
cat("Upregulated mRNAs (padj<", PADJ_THR, ", log2FC>=", LFC_UP, "): ",
    nrow(deg_up), "\n", sep = "")
cat("Top up-genes:", paste(head(deg_up$symbol, 15), collapse = ", "), "\n")

up_genes <- unique(deg_up$symbol)
cat("Unique up-gene symbols:", length(up_genes), "\n\n")

cat("Loading DEMs from", DEM_FILE, "...\n")
dem <- read.csv(DEM_FILE, stringsAsFactors = FALSE)
cat("Total downregulated miRNAs:", nrow(dem), "\n")

dem$miRNA_clean <- sub("\\|.*$", "", dem$miRNA)
expand_paired <- function(x) {
  if (grepl("\\+", x)) return(trimws(strsplit(x, "\\+")[[1]]))
  x
}
down_miRNAs <- unique(unlist(lapply(dem$miRNA_clean, expand_paired)))
cat("Unique downregulated miRNAs (after expanding paired probes):",
    length(down_miRNAs), "\n")
cat("List:\n"); print(down_miRNAs)

# ---- 2. Query multiMiR -----------------------------------------------------
cat("\nQuerying multiMiR for", length(down_miRNAs), "miRNAs against",
    length(up_genes), "target genes...\n\n")

mm <- suppressWarnings(get_multimir(
  org        = "hsa",
  mirna      = down_miRNAs,
  target     = up_genes,
  table      = "all",
  summary    = TRUE,
  use.tibble = FALSE
))

res_all <- mm@data
cat("Raw interaction rows returned:", nrow(res_all), "\n")

# Restrict to the 11 methodology databases and drop NA target symbols
res_all <- subset(res_all,
                  database %in% METHOD_DBS &
                  !is.na(target_symbol) & target_symbol != "")
cat("Rows after restricting to 11 methodology DBs and dropping NA targets:",
    nrow(res_all), "\n")
cat("Types present:\n"); print(table(res_all$type))
cat("Databases present:\n"); print(table(res_all$database))

# ---- 3. Separate validated vs predicted, summarize per pair ----------------
res_all$pair <- paste(res_all$mature_mirna_id, res_all$target_symbol,
                      sep = "::")

# Use multiMiR's own 'type' column, not database-list matching.
# type == "validated" -> experimentally validated targets.
# type == "predicted" -> sequence-based predictions.
# (type == "disease.drug" rows are miRNA-disease associations, not target
# validation, so they are excluded from both.)
res_val  <- subset(res_all, type == "validated")
res_pred <- subset(res_all, type == "predicted")

val_summary <- res_val %>%
  dplyr::group_by(pair, mature_mirna_id, target_symbol) %>%
  dplyr::summarise(
    n_val_dbs = dplyr::n_distinct(database),
    val_dbs   = paste(sort(unique(database)), collapse = ";"),
    val_methods = paste(sort(unique(experiment[!is.na(experiment) &
                                               experiment != ""])),
                        collapse = ";"),
    val_pubmed  = paste(sort(unique(pubmed_id[!is.na(pubmed_id) &
                                              pubmed_id != ""])),
                        collapse = ";"),
    .groups = "drop"
  )

pred_summary <- res_pred %>%
  dplyr::group_by(pair, mature_mirna_id, target_symbol) %>%
  dplyr::summarise(
    n_pred_dbs = dplyr::n_distinct(database),
    pred_dbs   = paste(sort(unique(database)), collapse = ";"),
    .groups = "drop"
  )

master <- dplyr::full_join(pred_summary, val_summary,
                           by = c("pair","mature_mirna_id","target_symbol"))
master$n_pred_dbs[is.na(master$n_pred_dbs)] <- 0
master$n_val_dbs[is.na(master$n_val_dbs)]   <- 0
master$pred_dbs[is.na(master$pred_dbs)]     <- ""
master$val_dbs[is.na(master$val_dbs)]       <- ""
master$validated      <- master$n_val_dbs > 0
master$high_confidence <- master$n_pred_dbs >= MIN_DBS_PRED | master$validated

# Merge in logFC / padj for ranking
master <- master %>%
  dplyr::left_join(
    deg_up[, c("symbol","log2FoldChange","padj")] %>%
      dplyr::rename(target_symbol = symbol,
                    mRNA_log2FC   = log2FoldChange,
                    mRNA_padj     = padj),
    by = "target_symbol") %>%
  dplyr::left_join(
    dem[, c("miRNA_clean","logFC","adj.P.Val")] %>%
      dplyr::rename(mature_mirna_id = miRNA_clean,
                    miRNA_logFC     = logFC,
                    miRNA_padj      = adj.P.Val),
    by = "mature_mirna_id")

master <- master[order(-master$validated, -master$n_pred_dbs,
                       -master$mRNA_log2FC), ]

cat("\nTotal unique miRNA-target pairs:", nrow(master), "\n")
cat("  Validated:", sum(master$validated), "\n")
cat("  High-confidence predicted (>=", MIN_DBS_PRED, " DBs): ",
    sum(master$n_pred_dbs >= MIN_DBS_PRED), "\n", sep = "")
cat("  High-confidence (validated OR >=", MIN_DBS_PRED, " predicted): ",
    sum(master$high_confidence), "\n", sep = "")

write.csv(master,
          file.path(OUT_DIR, "miRNA-target_all_pairs.csv"), row.names = FALSE)
write.csv(subset(master, validated),
          file.path(OUT_DIR, "miRNA-target_validated.csv"), row.names = FALSE)
write.csv(subset(master, high_confidence),
          file.path(OUT_DIR, "miRNA-target_high_confidence.csv"),
          row.names = FALSE)
write.csv(subset(master, n_pred_dbs > 0),
          file.path(OUT_DIR, "miRNA-target_predicted.csv"), row.names = FALSE)

# ---- 4. Per-miRNA and per-target summaries ---------------------------------
hc <- subset(master, high_confidence)

per_mirna <- hc %>%
  dplyr::group_by(mature_mirna_id) %>%
  dplyr::summarise(
    n_targets = dplyr::n_distinct(target_symbol),
    targets   = paste(sort(unique(target_symbol)), collapse = ", "),
    .groups = "drop"
  ) %>%
  dplyr::arrange(-n_targets)

per_target <- hc %>%
  dplyr::group_by(target_symbol) %>%
  dplyr::summarise(
    n_miRNAs = dplyr::n_distinct(mature_mirna_id),
    miRNAs   = paste(sort(unique(mature_mirna_id)), collapse = ", "),
    mRNA_log2FC = dplyr::first(mRNA_log2FC),
    .groups = "drop"
  ) %>%
  dplyr::arrange(-n_miRNAs)

cat("\n=== Per-miRNA summary (high-confidence targets) ===\n")
print(as.data.frame(per_mirna), row.names = FALSE)
cat("\n=== Per-target summary (high-confidence) ===\n")
print(as.data.frame(per_target), row.names = FALSE)

write.csv(per_mirna,  file.path(OUT_DIR, "summary_per_miRNA.csv"),
          row.names = FALSE)
write.csv(per_target, file.path(OUT_DIR, "summary_per_target.csv"),
          row.names = FALSE)

cat("\n=== DONE ===\n")
cat("Output folder:", normalizePath(OUT_DIR), "\n")
print(list.files(OUT_DIR))
