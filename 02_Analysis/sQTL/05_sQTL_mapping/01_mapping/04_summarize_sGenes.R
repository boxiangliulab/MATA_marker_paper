#!/usr/bin/env Rscript
# sGene tables.
#
# Summarise the significant sQTL of step 03 at the gene level: for each depot
# the genes carrying at least one significant sQTL ("sGenes"), with the
# number of sQTL, splice sites and cell types behind each one, plus the
# overlap of the sGenes between the two depots.

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
})

cat("=== sGene summarization ===\n\n")

# --- Paths (relative to this script's directory) -----------------------------
output_dir <- "../../data/05_sQTL_mapping/results"

sqtl_file       <- file.path(output_dir, "All_sQTL_SAT_and_VAT_with_beta.tsv")
output_sat      <- file.path(output_dir, "sGenes_SAT.tsv")
output_vat      <- file.path(output_dir, "sGenes_VAT.tsv")
output_combined <- file.path(output_dir, "sGenes_summary.tsv")

FDR_THRESHOLD <- 0.05

if (!file.exists(sqtl_file)) {
  cat(sprintf("ERROR: sQTL file not found: %s\n", sqtl_file))
  cat("Please run 03_summarize_sQTL.R first.\n")
  quit(status = 1)
}

# --- Read the combined sQTL table -------------------------------------------
cat("Reading sQTL summary file...\n")
sqtl_data <- fread(sqtl_file, header = TRUE)
cat(sprintf("Loaded %d sQTL associations\n", nrow(sqtl_data)))

# Filter for significant sQTLs (FDR < 0.05)
sqtl_sig <- sqtl_data %>%
  filter(FDR < FDR_THRESHOLD) %>%
  filter(!is.na(gene))

cat(sprintf("Significant sQTLs (FDR < %.2f): %d\n\n",
            FDR_THRESHOLD, nrow(sqtl_sig)))

# One row per gene with at least one significant sQTL in that depot
summarize_sgenes <- function(data, depot_name) {
  cat(sprintf("Processing %s...\n", depot_name))

  depot_data <- data %>%
    filter(depot == depot_name)

  if (nrow(depot_data) == 0) {
    cat(sprintf("  No significant sQTLs found for %s\n", depot_name))
    return(NULL)
  }

  sgenes <- depot_data %>%
    group_by(gene, gene_id, gene_type) %>%
    summarise(
      n_sQTLs = n(),
      n_splice_sites = n_distinct(splice_site),
      n_celltypes = n_distinct(celltype),
      celltypes = paste(sort(unique(celltype)), collapse = ";"),
      min_FDR = min(FDR, na.rm = TRUE),
      min_pvalue = min(p.value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(min_FDR, min_pvalue)

  cat(sprintf("  Unique sGenes: %d\n", nrow(sgenes)))
  cat(sprintf("  Total sQTLs: %d\n", sum(sgenes$n_sQTLs)))
  cat(sprintf("  Total splice sites: %d\n", sum(sgenes$n_splice_sites)))

  # Protein coding, lncRNA, ...
  gene_type_summary <- sgenes %>%
    group_by(gene_type) %>%
    summarise(n_genes = n(), .groups = "drop") %>%
    arrange(desc(n_genes))

  cat("\n  Gene type breakdown:\n")
  print(gene_type_summary)

  return(sgenes)
}

# --- SAT ---------------------------------------------------------------------
cat(paste(rep("=", 50), collapse = ""), "\n")
sgenes_sat <- summarize_sgenes(sqtl_sig, "SAT")

cat("\n")
cat(paste(rep("=", 50), collapse = ""), "\n")

# --- VAT ---------------------------------------------------------------------
sgenes_vat <- summarize_sgenes(sqtl_sig, "VAT")

# --- Save the per-depot tables -----------------------------------------------
if (!is.null(sgenes_sat)) {
  cat(sprintf("\nWriting SAT sGenes to: %s\n", output_sat))
  fwrite(sgenes_sat, output_sat, sep = "\t", quote = FALSE, row.names = FALSE)
}

if (!is.null(sgenes_vat)) {
  cat(sprintf("Writing VAT sGenes to: %s\n", output_vat))
  fwrite(sgenes_vat, output_vat, sep = "\t", quote = FALSE, row.names = FALSE)
}

# --- Counts per depot --------------------------------------------------------
cat("\n=== Combined Summary ===\n")

combined_summary <- data.frame(
  depot = c("SAT", "VAT"),
  n_sGenes = c(
    ifelse(!is.null(sgenes_sat), nrow(sgenes_sat), 0),
    ifelse(!is.null(sgenes_vat), nrow(sgenes_vat), 0)
  ),
  n_sQTLs = c(
    ifelse(!is.null(sgenes_sat), sum(sgenes_sat$n_sQTLs), 0),
    ifelse(!is.null(sgenes_vat), sum(sgenes_vat$n_sQTLs), 0)
  ),
  n_splice_sites = c(
    ifelse(!is.null(sgenes_sat), sum(sgenes_sat$n_splice_sites), 0),
    ifelse(!is.null(sgenes_vat), sum(sgenes_vat$n_splice_sites), 0)
  )
)

print(combined_summary)

cat(sprintf("\nWriting combined summary to: %s\n", output_combined))
fwrite(combined_summary, output_combined, sep = "\t", quote = FALSE, row.names = FALSE)

# --- Overlap of the sGenes between the depots --------------------------------
if (!is.null(sgenes_sat) && !is.null(sgenes_vat)) {
  sat_genes <- sgenes_sat$gene_id
  vat_genes <- sgenes_vat$gene_id

  overlap <- intersect(sat_genes, vat_genes)
  sat_specific <- setdiff(sat_genes, vat_genes)
  vat_specific <- setdiff(vat_genes, sat_genes)

  cat("\n=== Overlap Analysis ===\n")
  cat(sprintf("sGenes in both SAT and VAT: %d\n", length(overlap)))
  cat(sprintf("sGenes specific to SAT: %d\n", length(sat_specific)))
  cat(sprintf("sGenes specific to VAT: %d\n", length(vat_specific)))

  cat("\n=== Top 10 sGenes by number of sQTLs ===\n")

  cat("\nSAT:\n")
  top_sat <- sgenes_sat %>%
    arrange(desc(n_sQTLs)) %>%
    select(gene, gene_type, n_sQTLs, n_splice_sites, n_celltypes) %>%
    head(10)
  print(top_sat)

  cat("\nVAT:\n")
  top_vat <- sgenes_vat %>%
    arrange(desc(n_sQTLs)) %>%
    select(gene, gene_type, n_sQTLs, n_splice_sites, n_celltypes) %>%
    head(10)
  print(top_vat)
}

cat("\n=== DONE ===\n")
cat("Output files:\n")
cat(sprintf("  SAT sGenes: %s\n", output_sat))
cat(sprintf("  VAT sGenes: %s\n", output_vat))
cat(sprintf("  Summary: %s\n", output_combined))
