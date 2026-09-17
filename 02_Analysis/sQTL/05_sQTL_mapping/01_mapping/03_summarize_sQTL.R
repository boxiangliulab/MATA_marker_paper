#!/usr/bin/env Rscript
# Combined table of significant sQTL.
#
# Collect the significant sQTL of every cell type and depot into the single
# table that the downstream analyses use, adding the gene each splice site
# belongs to, the variant coordinates and the number of donors behind each
# cell type, and print the summary statistics reported in the manuscript.
#
# The FDR values are the ones computed in step 02; nothing is re-tested here.

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
})

cat("=== sQTL summarization (using the pre-calculated FDR) ===\n\n")

# --- Paths (relative to this script's directory) -----------------------------
results_dir    <- "../../data/05_sQTL_mapping/results"
pheno_dir      <- "../../data/04_phenotype"
# Splice site to gene annotation: gene_id, splice_site, gene_type, gene_name
site2gene_file <- "../../data/reference/site2gene.processed.tsv"

output_file <- file.path(results_dir, "All_sQTL_SAT_and_VAT_with_beta.tsv")

FDR_THRESHOLD <- 0.05

# --- Splice site to gene annotation ------------------------------------------
cat("Reading site2gene annotation...\n")
site2gene <- fread(site2gene_file, header = TRUE)
colnames(site2gene) <- c("gene_id", "splice_site", "gene_type", "gene_name")
cat(sprintf("Loaded %d splice sites mapped to genes\n", nrow(site2gene)))

# Number of donors behind a cell type: the meta-cell columns of its phenotype
# matrix are named <donor>:<index>, so the donors are the unique prefixes.
get_donor_count <- function(celltype, depot) {
  pheno_file <- file.path(pheno_dir, "phenotype_output", celltype, depot,
                          "phenotype_output.noFilter_rename.filtered")
  if (!file.exists(pheno_file)) return(NA_integer_)

  header <- tryCatch(readLines(pheno_file, n = 1), error = function(e) NULL)
  if (is.null(header)) return(NA_integer_)

  cols <- strsplit(header, " ")[[1]]
  donor_ids <- unique(sapply(strsplit(cols, ":"), function(x) x[1]))
  donor_ids <- donor_ids[grepl("^H\\d+$", donor_ids)]
  length(donor_ids)
}

# --- Read every significant-sQTL file ----------------------------------------
fdr_files <- list.files(results_dir, pattern = "\\.result\\.sQTL\\.fdr05$",
                        full.names = TRUE)
cat(sprintf("Found %d result.sQTL.fdr05 files\n\n", length(fdr_files)))

all_results <- list()

for (file_path in fdr_files) {
  filename <- basename(file_path)

  # File names are <celltype>.<depot>.p<p>_g<g>.result.sQTL.fdr05; the cell
  # type may itself contain dots, so it is everything before the depot field.
  parts <- strsplit(filename, "\\.")[[1]]

  if ("SAT" %in% parts) {
    depot <- "SAT"
    depot_idx <- which(parts == "SAT")
  } else if ("VAT" %in% parts) {
    depot <- "VAT"
    depot_idx <- which(parts == "VAT")
  } else {
    cat(sprintf("Warning: could not extract depot from %s, skipping\n", filename))
    next
  }

  celltype <- paste(parts[1:(depot_idx - 1)], collapse = ".")
  cat(sprintf("Processing: %s - %s\n", celltype, depot))

  n_donors <- get_donor_count(celltype, depot)
  cat(sprintf("  Donor count: %s\n",
              ifelse(is.na(n_donors), "NA", as.character(n_donors))))

  df <- tryCatch(fread(file_path, header = TRUE), error = function(e) {
    cat(sprintf("  Error reading %s: %s\n", filename, e$message))
    NULL
  })

  if (is.null(df) || nrow(df) == 0) {
    cat("  Skipping empty file\n")
    next
  }

  df$celltype <- celltype
  df$depot <- depot
  df$N <- n_donors

  cat(sprintf("  Total splice sites: %d, significant (FDR < %.2f): %d\n",
              nrow(df), FDR_THRESHOLD,
              sum(df$global_FDR < FDR_THRESHOLD, na.rm = TRUE)))

  all_results[[length(all_results) + 1]] <- df
}

# --- Combine ------------------------------------------------------------------
cat("\n=== Combining all results ===\n")
if (length(all_results) == 0) {
  cat("ERROR: no results were successfully processed!\n")
  quit(status = 1)
}

sqtl_combined <- rbindlist(all_results, fill = TRUE)
cat(sprintf("Total sQTLs: %d\n", nrow(sqtl_combined)))

# Variant IDs are chr:pos:ref:alt
cat("Parsing SNP information...\n")
sqtl_combined <- sqtl_combined %>%
  mutate(CHR     = sub("^(chr[^:]+):.*", "\\1", SNP),
         CHR_num = as.integer(sub("chr", "", CHR)),
         POS     = as.integer(sub("^chr[^:]+:(\\d+):.*", "\\1", SNP)),
         Allele1 = sub("^chr[^:]+:\\d+:([^:]+):.*", "\\1", SNP),
         Allele2 = sub("^chr[^:]+:\\d+:[^:]+:(.+)$", "\\1", SNP),
         MarkerID = SNP)

cat("Merging with site2gene annotation...\n")
sqtl_annotated <- sqtl_combined %>%
  left_join(site2gene %>% select(splice_site, gene_id, gene_type, gene_name),
            by = "splice_site")

missing_annot <- sum(is.na(sqtl_annotated$gene_name))
if (missing_annot > 0) {
  cat(sprintf("Warning: %d sQTLs have missing gene annotations\n", missing_annot))
}

cat("Formatting output table...\n")
sqtl_final <- sqtl_annotated %>%
  mutate(pos     = paste(CHR_num, POS, sep = ":"),
         gene    = gene_name,
         sQTL    = paste(splice_site, SNP, sep = "-"),
         FDR     = global_FDR,
         BETA    = beta,
         SE      = se,
         p.value = pvalue) %>%
  select(pos, CHR, POS, MarkerID, Allele1, Allele2,
         BETA, SE, p.value, N,
         depot, celltype, gene, FDR, sQTL,
         gene_id, gene_type, splice_site) %>%
  arrange(FDR, p.value)

cat(sprintf("\nWriting output to: %s\n", output_file))
fwrite(sqtl_final, output_file, sep = "\t", quote = FALSE, row.names = FALSE,
       na = "NA")

# --- Summary statistics -------------------------------------------------------
cat("\n=== Summary statistics ===\n")
cat(sprintf("Total sQTLs: %d\n", nrow(sqtl_final)))
cat(sprintf("Unique splice sites: %d\n", length(unique(sqtl_final$splice_site))))
cat(sprintf("Unique genes: %d\n",
            length(unique(sqtl_final$gene[!is.na(sqtl_final$gene)]))))
cat(sprintf("Unique SNPs: %d\n", length(unique(sqtl_final$MarkerID))))

cat("\nBreakdown by depot:\n")
print(sqtl_final %>%
        group_by(depot) %>%
        summarise(total = n(),
                  significant = sum(FDR < FDR_THRESHOLD, na.rm = TRUE),
                  pct_sig = 100 * significant / total,
                  .groups = "drop"))

cat("\nTop 15 cell types by significant sQTLs:\n")
print(sqtl_final %>%
        filter(FDR < FDR_THRESHOLD) %>%
        group_by(celltype, depot) %>%
        summarise(n_sqtl = n(), .groups = "drop") %>%
        arrange(desc(n_sqtl)) %>%
        head(15))

cat("\nTop 10 genes by number of significant sQTLs:\n")
print(sqtl_final %>%
        filter(FDR < FDR_THRESHOLD, !is.na(gene)) %>%
        group_by(gene) %>%
        summarise(n_sqtl = n(), .groups = "drop") %>%
        arrange(desc(n_sqtl)) %>%
        head(10))

cat("\nSample size (N) by cell type:\n")
print(sqtl_final %>%
        group_by(celltype, depot) %>%
        summarise(N = first(N), .groups = "drop") %>%
        arrange(depot, celltype))

cat("\n=== DONE ===\n")
cat(sprintf("Output saved to: %s\n", output_file))
