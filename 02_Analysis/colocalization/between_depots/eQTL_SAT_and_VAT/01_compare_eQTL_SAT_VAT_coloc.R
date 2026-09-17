#!/usr/bin/env Rscript
# eQTL colocalization SAT vs VAT.
#
# The eQTL counterpart of 01_sQTL_SAT_and_VAT: per cell type, decide
# whether the cis-eQTL signals seen in both depots are driven by the same
# causal variant. Genes with a candidate signal (p < 1e-3) in both SAT and
# VAT are tested with coloc.abf; genes with PP.H4 > 0.8 are called shared,
# the remaining candidate genes are reported as depot-specific.
#
# coloc is run on p-values and MAF (no sdY), with the per-depot sample sizes
# taken from the SAIGE-QTL output.

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(readr)
  library(purrr)
  library(tidyr)
  library(coloc)
  library(data.table)
})

# ====== Configuration ======
# Paths are relative to this script's directory.
# SAIGE-QTL eQTL results: candidate loci and the per-gene cis results
eqtl_loci_dir    <- "../../data/eQTL/loci"
eqtl_results_dir <- "../../data/eQTL/cis_results"
# Cohort MAF per population
maf_file         <- "../../data/reference/MAF_3pop_merged.tsv"
out_dir          <- "../../data/07_colocalization/eQTL_SAT_and_VAT"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Colocalization thresholds
PP4_THRESHOLD <- 0.8   # posterior probability of a shared causal variant
MIN_VARIANTS  <- 50    # minimum number of variants shared by the two depots

# ====== Load MAF data ======
cat("Loading MAF data...\n")
MAF <- fread(maf_file)
MAF$variant_id <- gsub("^chr", "", MAF$variant_id)
MAF[, MAF_Chinese := Chinese]
MAF[, MAF_Indian  := Indian]
MAF[, MAF_Malay   := Malay]
MAF[, MAF_mean := rowMeans(.SD, na.rm = TRUE),
    .SDcols = c("MAF_Chinese", "MAF_Indian", "MAF_Malay")]
MAF <- MAF[, .(variant_id, rsID, variant, MAF_Chinese, MAF_Indian, MAF_Malay,
               MAF_mean)]

# ====== Helper functions ======

# Read the candidate eQTL loci of one depot (genes with p < 1e-3).
# The file is raw SAIGE-QTL output without a header.
read_eqtl_loci <- function(file_path, depot) {
  if (!file.exists(file_path)) return(NULL)

  df <- fread(file_path, header = FALSE)
  colnames(df) <- c("CHR", "POS", "variant_id", "Allele1", "Allele2",
                    "AC_Allele2", "AF_Allele2", "MissingRate", "BETA", "SE",
                    "Tstat", "var", "p.value", "p.value.NA", "Is.SPA", "N",
                    "gene")

  # Variant IDs are reduced to chr:pos (no "chr" prefix) so that they match
  # the MAF table and the two depots can be merged.
  df$variant_id <- gsub("^chr", "", df$variant_id)
  df$variant_id <- sub("(:[^:]+){2}$", "", df$variant_id)
  df$depot <- depot

  df
}

# Cache of the per-gene result files (each gene is read once per depot)
gene_cache <- new.env()

# Read all cis variants tested for one gene in one depot
read_all_variants_for_gene <- function(cell_type, depot, gene, eqtl_results_dir) {
  cache_key <- paste0(cell_type, ".", depot, ".", gene)

  if (!exists(cache_key, envir = gene_cache)) {
    gene_file <- file.path(eqtl_results_dir, depot, cell_type,
                           paste0(gene, "_cis"))

    if (!file.exists(gene_file)) {
      warning("Gene file not found: ", gene_file)
      return(NULL)
    }

    cat("    Loading full results for", gene, "in", depot, "...\n")

    df <- fread(gene_file, sep = "\t", header = TRUE)
    df$gene <- gene
    df <- df[, .(gene, MarkerID, CHR, POS, BETA, SE, p.value, N)]

    # chr:pos variant IDs, as for the loci files above
    df[, variant_id := gsub("^chr", "", MarkerID)]
    df[, variant_id := sapply(strsplit(variant_id, ":"),
                              function(x) paste(x[1], x[2], sep = ":"))]

    assign(cache_key, df, envir = gene_cache)
  }

  get(cache_key, envir = gene_cache)
}

# Colocalize the SAT and VAT signals of one gene
run_coloc_for_gene <- function(gene_id, cell_type, eqtl_results_dir,
                               progress_str = "") {

  cat(progress_str, "Coloc:", gene_id, "\n")

  sat_full <- read_all_variants_for_gene(cell_type, "SAT", gene_id,
                                         eqtl_results_dir)
  vat_full <- read_all_variants_for_gene(cell_type, "VAT", gene_id,
                                         eqtl_results_dir)

  if (is.null(sat_full) || is.null(vat_full)) {
    return(list(PP.H4.abf = NA, nsnps = 0, status = "missing_data"))
  }
  if (nrow(sat_full) == 0 || nrow(vat_full) == 0) {
    return(list(PP.H4.abf = NA, nsnps = 0, status = "no_variants"))
  }

  # Variants tested in both depots
  merged <- inner_join(
    sat_full %>% select(variant_id, BETA, SE, p.value, N),
    vat_full %>% select(variant_id, BETA, SE, p.value, N),
    by = "variant_id",
    suffix = c("_SAT", "_VAT")
  )

  if (nrow(merged) < MIN_VARIANTS) {
    return(list(PP.H4.abf = NA, nsnps = nrow(merged),
                status = "insufficient_variants"))
  }

  # Attach the cohort MAF
  merged_dt <- as.data.table(merged)
  merged_dt <- merge(merged_dt, MAF, by = "variant_id", all.x = TRUE)
  merged_dt <- merged_dt[!is.na(MAF_Chinese) & MAF_Chinese > 0 & MAF_Chinese < 1]

  if (nrow(merged_dt) < MIN_VARIANTS) {
    return(list(PP.H4.abf = NA, nsnps = nrow(merged_dt),
                status = "insufficient_variants_after_maf"))
  }

  # Drop incomplete records
  merged_dt <- merged_dt[!is.na(p.value_SAT) & !is.na(p.value_VAT) &
                         !is.na(BETA_SAT)    & !is.na(BETA_VAT) &
                         !is.na(SE_SAT)      & !is.na(SE_VAT) &
                         SE_SAT > 0 & SE_VAT > 0]

  # One record per variant; keep the most significant when duplicated
  merged_dt <- merged_dt[order(p.value_SAT)][!duplicated(variant_id)]

  if (nrow(merged_dt) < MIN_VARIANTS) {
    return(list(PP.H4.abf = NA, nsnps = nrow(merged_dt),
                status = "insufficient_variants_after_filter"))
  }

  # Sample sizes reported by SAIGE-QTL for this gene
  N_SAT <- unique(merged_dt$N_SAT)[1]
  N_VAT <- unique(merged_dt$N_VAT)[1]

  cat("    Using", nrow(merged_dt), "variants, N_SAT =", N_SAT,
      ", N_VAT =", N_VAT, "\n")

  # coloc datasets: p-values + MAF (no sdY), quantitative traits
  dataset1 <- list(pvalues = merged_dt$p.value_SAT, N = N_SAT,
                   type = "quant", snp = merged_dt$variant_id)
  dataset2 <- list(pvalues = merged_dt$p.value_VAT, N = N_VAT,
                   type = "quant", snp = merged_dt$variant_id)

  tryCatch({
    result <- coloc.abf(dataset1, dataset2, MAF = merged_dt$MAF_Chinese)
    list(PP.H0.abf = result$summary["PP.H0.abf"],
         PP.H1.abf = result$summary["PP.H1.abf"],
         PP.H2.abf = result$summary["PP.H2.abf"],
         PP.H3.abf = result$summary["PP.H3.abf"],
         PP.H4.abf = result$summary["PP.H4.abf"],
         nsnps = result$summary["nsnps"],
         status = "success")
  }, error = function(e) {
    cat("    ERROR:", e$message, "\n")
    list(PP.H4.abf = NA, nsnps = nrow(merged_dt), status = "coloc_error")
  })
}

# ====== Analysis of one cell type ======

analyse_cell_type <- function(CELL_TYPE) {


  cat("\n========================================\n")
  cat("Cell type:", CELL_TYPE, "\n")
  cat("========================================\n")

  sat_file <- file.path(eqtl_loci_dir, "SAT", "1e-3",
                        paste0("SAT.", CELL_TYPE, ".1e-3.loci"))
  vat_file <- file.path(eqtl_loci_dir, "VAT", "1e-3",
                        paste0("VAT.", CELL_TYPE, ".1e-3.loci"))

  if (!file.exists(sat_file)) stop("SAT file not found: ", sat_file)
  if (!file.exists(vat_file)) stop("VAT file not found: ", vat_file)

  sat_eqtl <- read_eqtl_loci(sat_file, "SAT")
  vat_eqtl <- read_eqtl_loci(vat_file, "VAT")

  if (is.null(sat_eqtl) || is.null(vat_eqtl)) stop("Failed to read eQTL data")

  sat_genes <- unique(sat_eqtl$gene)
  vat_genes <- unique(vat_eqtl$gene)

  cat("  SAT genes:", length(sat_genes), "\n")
  cat("  VAT genes:", length(vat_genes), "\n")

  # Genes with a candidate signal in both depots
  shared_genes <- intersect(sat_genes, vat_genes)
  cat("  Shared genes (candidates):", length(shared_genes), "\n")

  ct_out <- file.path(out_dir, CELL_TYPE)
  dir.create(ct_out, recursive = TRUE, showWarnings = FALSE)

  if (length(shared_genes) == 0) {
    cat("  No shared genes, saving empty results\n")

    write_tsv(tibble(gene = character(0)), file.path(ct_out, "colocalized.tsv"))
    write_tsv(tibble(gene = sat_genes), file.path(ct_out, "SAT_specific.tsv"))
    write_tsv(tibble(gene = vat_genes), file.path(ct_out, "VAT_specific.tsv"))

    summary_df <- tibble(cell_type = CELL_TYPE,
                         n_SAT = length(sat_genes),
                         n_VAT = length(vat_genes),
                         n_shared_candidates = 0,
                         n_colocalized = 0,
                         n_SAT_specific = length(sat_genes),
                         n_VAT_specific = length(vat_genes))
    write_tsv(summary_df, file.path(ct_out, "summary.tsv"))

    return(summary_df)
  }

  cat("  Running colocalization...\n")
  n_shared <- length(shared_genes)

  coloc_results <- map_dfr(seq_along(shared_genes), function(i) {
    gene_id <- shared_genes[i]
    progress <- sprintf("    [%d/%d] ", i, n_shared)
    res <- run_coloc_for_gene(gene_id, CELL_TYPE, eqtl_results_dir, progress)

    tibble(gene = gene_id,
           PP.H0.abf = res$PP.H0.abf %||% NA,
           PP.H1.abf = res$PP.H1.abf %||% NA,
           PP.H2.abf = res$PP.H2.abf %||% NA,
           PP.H3.abf = res$PP.H3.abf %||% NA,
           PP.H4.abf = res$PP.H4.abf %||% NA,
           nsnps = res$nsnps,
           status = res$status)
  })

  # ====== Classify the genes ======
  colocalized_genes <- coloc_results %>%
    filter(!is.na(PP.H4.abf), PP.H4.abf > PP4_THRESHOLD) %>%
    pull(gene)

  not_colocalized <- setdiff(shared_genes, colocalized_genes)
  sat_specific <- setdiff(sat_genes, c(colocalized_genes, not_colocalized))
  vat_specific <- setdiff(vat_genes, c(colocalized_genes, not_colocalized))

  cat("  Colocalized (PP.H4 >", PP4_THRESHOLD, "):", length(colocalized_genes), "\n")
  cat("  SAT-specific:", length(sat_specific), "\n")
  cat("  VAT-specific:", length(vat_specific), "\n")

  write_tsv(tibble(gene = colocalized_genes), file.path(ct_out, "colocalized.tsv"))
  write_tsv(tibble(gene = sat_specific), file.path(ct_out, "SAT_specific.tsv"))
  write_tsv(tibble(gene = vat_specific), file.path(ct_out, "VAT_specific.tsv"))
  write_tsv(coloc_results, file.path(ct_out, "coloc_results_full.tsv"))

  summary_df <- tibble(cell_type = CELL_TYPE,
                       n_SAT = length(sat_genes),
                       n_VAT = length(vat_genes),
                       n_shared_candidates = length(shared_genes),
                       n_colocalized = length(colocalized_genes),
                       n_SAT_specific = length(sat_specific),
                       n_VAT_specific = length(vat_specific))
  write_tsv(summary_df, file.path(ct_out, "summary.tsv"))

  cat("\nResults saved to:", ct_out, "\n")
  cat("========================================\n")

  return(summary_df)
}

# ====== Entry point ======
# With one or more cell types on the command line only those are processed;
# with no argument every cell type tested in both depots is processed in turn
# and the per-cell-type summaries are combined.

ref_dir <- "../../reference"

celltypes_both_depots <- function() {
  read_list <- function(depot) {
    f <- file.path(ref_dir, paste0("celltypes.QTL.", depot, ".txt"))
    if (!file.exists(f)) stop("Cell-type list not found: ", f)
    read.table(f, sep = "\t", stringsAsFactors = FALSE)[[1]]
  }
  intersect(read_list("SAT"), read_list("VAT"))
}

args <- commandArgs(trailingOnly = TRUE)
run_all <- length(args) == 0
cell_types <- if (run_all) celltypes_both_depots() else args

summaries <- lapply(cell_types, analyse_cell_type)

if (run_all) {
  summary_all <- bind_rows(summaries)
  summary_file <- file.path(out_dir, "summary_all_celltypes.tsv")
  write_tsv(summary_all, summary_file)
  cat("\nCombined summary written to:", summary_file, "\n")
}
