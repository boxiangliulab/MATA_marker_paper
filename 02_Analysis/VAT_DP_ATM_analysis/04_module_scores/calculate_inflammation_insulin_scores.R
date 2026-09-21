#!/usr/bin/env Rscript
# Purpose: Calculate inflammatory and insulin-signalling module scores by cell
# type and summarize them at donor level.
# Inputs: SECTION4_METADATA and raw count matrices under SECTION4_COUNT_ROOT.
# Outputs: Per-cell score RDS files and donor-mean TSV files under
# SECTION4_OUTPUT/module_scores.
# Method: Seurat AddModuleScore using the gene sets defined below; donors with
# fewer than 10 cells of a cell type are excluded from donor-level summaries.
# Usage: Rscript calculate_inflammation_insulin_scores.R

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(msigdbr)
  library(Seurat)
})

set.seed(123)
options(error = function() { traceback(20); q(status = 1) })

# =========================
# 0) Paths & settings
# =========================
metadata_path  <- Sys.getenv("SECTION4_METADATA")
count_base     <- Sys.getenv("SECTION4_COUNT_ROOT")

out_root <- file.path(Sys.getenv("SECTION4_OUTPUT"), "module_scores")
out_rds  <- file.path(out_root, "rds")   # per-cell raw scores (RDS)
out_tsv  <- file.path(out_root, "tsv")   # donor-mean per celltype (TSV)

dir.create(out_rds, recursive = TRUE, showWarnings = FALSE)
dir.create(out_tsv, recursive = TRUE, showWarnings = FALSE)

donor_col     <- "donor_id"
cellstate_col <- "Annotation_level2"
min_cells     <- 10

# =========================
# 1. Load metadata and enumerate level-2 cell annotations.
# =========================
metadata  <- read.table(metadata_path, sep = "\t", header = TRUE, check.names = FALSE)

# Derive cell types directly from the frozen metadata.
celltypes <- metadata %>%
  pull(Annotation_level2) %>%
  unique() %>%
  na.omit() %>%
  as.character()

cat("Total celltypes:", length(celltypes), "\n")
cat(celltypes, sep = "\n")

# =========================
# 2) Define gene sets
# =========================
gs_H  <- msigdbr(species = "Homo sapiens", collection = "H")
gs_C2 <- msigdbr(species = "Homo sapiens", collection = "C2")
gs_C5 <- msigdbr(species = "Homo sapiens", collection = "C5")

get_genes <- function(gs_df, set_name) {
  gs_df %>% filter(gs_name == set_name) %>% pull(gene_symbol) %>% unique()
}

all_sets <- list(
  # Insulin-signalling gene sets.
  MS_KEGG_INSULIN        = get_genes(gs_C2, "KEGG_INSULIN_SIGNALING_PATHWAY"),
  MS_KEGG_T2DM           = get_genes(gs_C2, "KEGG_TYPE_II_DIABETES_MELLITUS"),
  MS_REACTOME_IR         = get_genes(gs_C2, "REACTOME_SIGNALING_BY_INSULIN_RECEPTOR"),
  MS_REACTOME_IR_CASCADE = get_genes(gs_C2, "REACTOME_INSULIN_RECEPTOR_SIGNALLING_CASCADE"),
  MS_GOBP_IR             = get_genes(gs_C5, "GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY"),
  MS_GOBP_INSULIN        = get_genes(gs_C5, "GOBP_RESPONSE_TO_INSULIN"),
  # Inflammatory gene sets.
  MS_TNFA_NFKB           = get_genes(gs_H,  "HALLMARK_TNFA_SIGNALING_VIA_NFKB"),
  MS_INFLAM              = get_genes(gs_H,  "HALLMARK_INFLAMMATORY_RESPONSE"),
  MS_IL6_JAK             = get_genes(gs_H,  "HALLMARK_IL6_JAK_STAT3_SIGNALING"),
  MS_IFNG                = get_genes(gs_H,  "HALLMARK_INTERFERON_GAMMA_RESPONSE"),
  MS_IL1B                = get_genes(gs_C5, "GOBP_INTERLEUKIN_1_BETA_PRODUCTION"),
  MS_IL10                = get_genes(gs_C5, "GOBP_INTERLEUKIN_10_PRODUCTION")
)

cat("Gene set sizes:\n")
sapply(all_sets, length) %>% print()

score_cols <- names(all_sets)

# =========================
# 3) Helper: load count matrix → Seurat object (VAT only)
# =========================
generate_obj <- function(celltype, metadata) {
  path <- file.path(count_base, "VAT", "01.raw_count",
                    paste0(celltype, ".raw_count.tsv"))
  obj <- fread(path, sep = "\t", header = TRUE, check.names = FALSE)
  obj <- as.data.frame(obj)
  rownames(obj) <- obj$V1
  obj$V1 <- NULL
  obj <- t(obj) %>% as.matrix()

  seu <- CreateSeuratObject(obj)
  meta_sub <- metadata[metadata$barcode %in% rownames(seu@meta.data), , drop = FALSE]
  rownames(meta_sub) <- meta_sub$barcode
  seu <- AddMetaData(seu, meta_sub)
  return(seu)
}

# =========================
# 4) Main loop: per celltype
# =========================
for (celltype in celltypes) {
  cat("\n====", celltype, "====\n")

  # --- load VAT ---
  seu <- tryCatch(
    generate_obj(celltype, metadata),
    error = function(e) { message("VAT load failed: ", e$message); NULL }
  )
  if (is.null(seu)) { message("Skipping ", celltype); next }

  seu <- JoinLayers(seu, layers = "counts")
  seu <- NormalizeData(seu, assay = "RNA", verbose = FALSE)

  genes_in_obj <- rownames(seu)

  # --- compute module scores ---
  for (nm in score_cols) {
    gs_present <- intersect(all_sets[[nm]], genes_in_obj)
    cat(sprintf("  %s: %d / %d genes present\n", nm, length(gs_present), length(all_sets[[nm]])))

    if (length(gs_present) < 5) {
      seu[[nm]] <- NA_real_
      next
    }

    seu <- AddModuleScore(
      seu,
      features = list(gs_present),
      name     = paste0(nm, "_tmp_"),
      seed     = 123
    )
    seu[[nm]] <- seu[[paste0(nm, "_tmp_1")]]
    seu[[paste0(nm, "_tmp_1")]] <- NULL
  }

  # --- save per-cell raw scores as RDS ---
  cell_meta <- seu[[]] %>%
    tibble::rownames_to_column("cell_id") %>%
    select(cell_id,
           donor             = all_of(donor_col),
           Annotation_level2 = all_of(cellstate_col),
           all_of(score_cols))

  saveRDS(
    cell_meta,
    file = file.path(out_rds, paste0(celltype, "_module_scores.rds"))
  )
  cat("  Saved RDS:", celltype, "\n")

  # --- aggregate to donor level ---
  donor_mean <- cell_meta %>%
    group_by(donor, Annotation_level2) %>%
    summarise(
      n_cells = n(),
      across(
        all_of(score_cols),
        ~ ifelse(n() >= min_cells, mean(.x, na.rm = TRUE), NA_real_)
      ),
      .groups = "drop"
    ) %>%
    filter(Annotation_level2 == celltype)

  write.table(
    donor_mean,
    file      = file.path(out_tsv, paste0(celltype, "_donor_mean_MS.tsv")),
    sep       = "\t",
    quote     = FALSE,
    row.names = FALSE
  )
  cat("  Saved TSV:", celltype, "\n")

  rm(seu); gc()
}

cat("\n\nDONE. Outputs in:", out_root, "\n")
cat("  RDS (per-cell):   ", out_rds, "\n")
cat("  TSV (donor-mean): ", out_tsv, "\n")
