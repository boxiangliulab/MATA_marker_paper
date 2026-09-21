#!/usr/bin/env Rscript
# Purpose: Generate donor-level pseudobulk expression matrices by depot and
# cell type for the Section 4 macrophage and mesothelial analyses.
#
# Inputs:
#   1. Cell-level metadata TSV containing barcode, donor_id, and Depot.
#   2. SAT raw-count directory containing <cell_type>.raw_count.tsv files.
#   3. VAT raw-count directory containing <cell_type>.raw_count.tsv files.
#
# Outputs:
#   - Raw donor-level pseudobulk count matrices.
#   - TMM-normalized log2(CPM + 1) matrices.
#
# Method:
#   Counts are summed across nuclei from the same donor, depot, and cell type.
#   Donors must contribute at least 10 nuclei for the cell type. Genes retained
#   in the normalized matrix must be detected in at least 10% of donors.
#
# Usage:
#   Rscript generate_donor_pseudobulk.R METADATA SAT_COUNTS VAT_COUNTS OUT_DIR \
#     [DEPOTS] [CELL_TYPES]
#   DEPOTS and CELL_TYPES are optional comma-separated lists.
#   With no arguments, paths are read from SECTION4_METADATA,
#   SECTION4_SAT_COUNTS, SECTION4_VAT_COUNTS, and SECTION4_OUTPUT.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(Seurat)
  library(edgeR)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  meta_path <- Sys.getenv("SECTION4_METADATA")
  sat_dir   <- Sys.getenv("SECTION4_SAT_COUNTS")
  vat_dir   <- Sys.getenv("SECTION4_VAT_COUNTS")
  out_root  <- file.path(Sys.getenv("SECTION4_OUTPUT"), "pseudobulk")
} else if (length(args) < 4) {
  stop("Usage: Rscript generate_donor_pseudobulk.R METADATA SAT_COUNTS VAT_COUNTS OUT_DIR [DEPOTS] [CELL_TYPES]")
} else {
  meta_path <- args[[1]]
  sat_dir   <- args[[2]]
  vat_dir   <- args[[3]]
  out_root  <- args[[4]]
}

DEPOTS_TO_RUN <- if (length(args) >= 5) strsplit(args[[5]], ",", fixed = TRUE)[[1]] else "VAT"
CELLTYPES_TO_RUN <- if (length(args) >= 6) {
  strsplit(args[[6]], ",", fixed = TRUE)[[1]]
} else {
  c("TIM4+_ATM", "TIM4+_CD11c+_ATM", "Mesothelial", "LAM", "PVM")
}

MIN_FRAC_DONORS <- 0.10   # gene kept if expressed (>0) in >= this fraction of donors (within depot)

DONOR_MIN_NUCLEI <- 10    # per-depot, per-celltype donor filter
MIN_DONORS       <- 3     # after donor filter
MIN_TOTAL_NUCLEI <- 500   # after donor filter

raw_path <- function(depot, ct) file.path(
  if (depot == "SAT") sat_dir else vat_dir,
  paste0(ct, ".raw_count.tsv")
)

out_sat  <- file.path(out_root, "SAT")
out_vat  <- file.path(out_root, "VAT")
dir.create(out_sat, recursive = TRUE, showWarnings = FALSE)
dir.create(out_vat, recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# load metadata
# ----------------------------
metadata <- fread(meta_path, sep = "\t", header = TRUE) %>% as.data.frame()
need_cols <- c("barcode", "donor_id", "Depot")
miss <- setdiff(need_cols, colnames(metadata))
if (length(miss) > 0) stop("Missing metadata columns: ", paste(miss, collapse = ", "))

# ----------------------------
# helper: read raw_count.tsv -> Seurat
# ----------------------------
read_raw_to_seurat <- function(tsv_path, metadata) {

  obj <- fread(tsv_path, sep = "\t", header = TRUE, check.names = FALSE)
  obj <- as.data.frame(obj)

  if (!"V1" %in% colnames(obj)) stop("[ERROR] raw_count has no V1 column: ", tsv_path)

  rownames(obj) <- obj$V1
  obj$V1 <- NULL

  make_seu_try <- function(mat) {
    seu <- CreateSeuratObject(mat)
    cells <- colnames(seu)
    meta_sub <- metadata[match(cells, metadata$barcode), , drop = FALSE]
    list(seu = seu, meta_sub = meta_sub)
  }

  # try 1) assume cells x genes -> transpose
  mat1 <- t(as.matrix(obj))
  out1 <- make_seu_try(mat1)
  miss1 <- sum(is.na(out1$meta_sub$barcode))

  if (miss1 == 0) {
    rownames(out1$meta_sub) <- out1$meta_sub$barcode
    out1$meta_sub <- out1$meta_sub[colnames(out1$seu), , drop = FALSE]
    out1$seu <- AddMetaData(out1$seu, out1$meta_sub)
    return(out1$seu)
  }

  # try 2) assume genes x cells -> no transpose
  mat2 <- as.matrix(obj)
  out2 <- make_seu_try(mat2)
  miss2 <- sum(is.na(out2$meta_sub$barcode))

  if (miss2 == 0) {
    rownames(out2$meta_sub) <- out2$meta_sub$barcode
    out2$meta_sub <- out2$meta_sub[colnames(out2$seu), , drop = FALSE]
    out2$seu <- AddMetaData(out2$seu, out2$meta_sub)
    return(out2$seu)
  }

  cells1 <- colnames(out1$seu)
  ex1 <- head(cells1[is.na(out1$meta_sub$barcode)], 5)

  cells2 <- colnames(out2$seu)
  ex2 <- head(cells2[is.na(out2$meta_sub$barcode)], 5)

  stop(
    "[ERROR] Cell barcodes could not be matched to the metadata.\n",
    " - transpose attempt: ", miss1, " unmatched (examples: ", paste(ex1, collapse=", "), ")\n",
    " - untransposed attempt: ", miss2, " unmatched (examples: ", paste(ex2, collapse=", "), ")\n",
    "Check the matrix orientation and barcode format.\n",
    "file: ", tsv_path
  )
}

# ----------------------------
# normalize (TMM -> CPM -> log2(CPM+1))
# note: applies expressed-in-donors filter before normalization
# ----------------------------
log2cpm_tmm <- function(mat_counts, min_frac_donors = 0.10) {

  mat_counts <- as.matrix(mat_counts)

  if (ncol(mat_counts) >= 1) {
    keep_genes <- rowSums(mat_counts > 0) >= (min_frac_donors * ncol(mat_counts))
    mat_counts <- mat_counts[keep_genes, , drop = FALSE]
  }

  dge <- DGEList(counts = mat_counts)
  dge <- calcNormFactors(dge, method = "TMM")

  norm_cpm <- cpm(dge, normalized.lib.sizes = TRUE, log = FALSE)
  log2cpm  <- log2(norm_cpm + 1)

  return(log2cpm)
}

# ----------------------------
# run one depot + one celltype
# ----------------------------
run_one_depot <- function(depot, celltype) {
  ct_out <- gsub("\\+", "plus", celltype)
  in_tsv <- raw_path(depot, celltype)
  if (!file.exists(in_tsv)) {
    message(paste0("[SKIP] ", depot, " missing raw_count for ", celltype))
    return(invisible(NULL))
  }

  message(paste0("[RUN] ", depot, " ", celltype))
  seu <- read_raw_to_seurat(in_tsv, metadata)

  # subset to correct depot (safety)
  seu <- subset(seu, Depot == depot)

  # donor nuclei >= DONOR_MIN_NUCLEI (within depot)
  df_donors <- as.data.frame(table(seu$donor_id))
  donors_keep <- as.character(df_donors$Var1[df_donors$Freq >= DONOR_MIN_NUCLEI])

  if (length(donors_keep) < MIN_DONORS) {
    message(paste0("[SKIP] ", depot, " ", celltype,
                   " donors < ", MIN_DONORS, " after nuclei>=", DONOR_MIN_NUCLEI))
    return(invisible(NULL))
  }

  seu <- subset(seu, donor_id %in% donors_keep)

  # total nuclei >= MIN_TOTAL_NUCLEI
  if (ncol(seu) < MIN_TOTAL_NUCLEI) {
    message(paste0("[SKIP] ", depot, " ", celltype,
                   " nuclei < ", MIN_TOTAL_NUCLEI, " after donor filter"))
    return(invisible(NULL))
  }

  # pseudobulk donor_id (raw counts, SUM across cells)
  mat_counts <- AggregateExpression(
    seu,
    assays   = "RNA",
    group.by = "donor_id",
    slot     = "counts"
  )$RNA

  out_dir <- if (depot == "SAT") out_sat else out_vat

  # Save unfiltered raw pseudobulk counts.
  out_counts <- file.path(out_dir, paste0(ct_out, "__", depot, "__donor_pseudobulk_COUNTS.tsv"))
  df_counts <- cbind(Gene = rownames(mat_counts), as.data.frame(mat_counts))
  fwrite(df_counts, out_counts, sep = "\t")

  # 2) save normalized log2CPM(TMM) (gene filter applied inside log2cpm_tmm)
  log2cpm <- log2cpm_tmm(mat_counts, min_frac_donors = MIN_FRAC_DONORS)
  out_log2cpm <- file.path(out_dir, paste0(ct_out, "__", depot, "__donor_pseudobulk_log2CPM_TMM.tsv"))
  df_log2cpm <- cbind(Gene = rownames(log2cpm), as.data.frame(log2cpm))
  fwrite(df_log2cpm, out_log2cpm, sep = "\t")

  message(paste0("[OK] saved COUNTS: ", out_counts))
  message(paste0("[OK] saved log2CPM(TMM): ", out_log2cpm))
  message(paste0("[INFO] ", depot, " ", celltype,
                 " donors=", ncol(mat_counts),
                 " genes(counts)=", nrow(mat_counts),
                 " genes(log2cpm)=", nrow(log2cpm)))

  invisible(NULL)
}

# ----------------------------
# main
# ----------------------------
DEPOTS_TO_RUN <- toupper(DEPOTS_TO_RUN)
if (!all(DEPOTS_TO_RUN %in% c("SAT", "VAT"))) stop("[ERROR] DEPOTS_TO_RUN must be SAT and/or VAT")

message(paste0("[INFO] depots    = ", paste(DEPOTS_TO_RUN, collapse = ", ")))
message(paste0("[INFO] celltypes = ", paste(CELLTYPES_TO_RUN, collapse = ", ")))

for (ct in CELLTYPES_TO_RUN) {
  message(paste0("========== ", ct, " =========="))
  for (dp in DEPOTS_TO_RUN) {
    run_one_depot(dp, ct)
  }
}

message("[DONE] selected celltypes finished.")
