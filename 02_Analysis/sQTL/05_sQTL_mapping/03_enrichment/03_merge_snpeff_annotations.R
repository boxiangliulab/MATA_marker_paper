#!/usr/bin/env Rscript
# Join the annotations to the sQTL.
#
# Parse the SnpEff-annotated VCF of step 02 and join the annotations back
# onto the significant sQTL, producing the annotated table that the
# enrichment test (step 06) reads. A variant with several SnpEff annotations
# contributes one row per annotation.

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
})

source("enrichment_functions.R")

# --- Paths (relative to this script's directory) -----------------------------
results_dir <- "../../data/05_sQTL_mapping/results"
snpeff_dir  <- "../../data/05_sQTL_mapping/SNPeff"

output_dir    <- snpeff_dir
annotated_vcf <- file.path(output_dir, "sQTL_variants_FDR05_annotated.vcf")
sqtl_file     <- file.path(results_dir,
                           "All_sQTL_SAT_and_VAT_with_beta.tsv")

FDR_THRESHOLD <- 0.05

cat("=== Joining SnpEff annotations to the significant sQTL ===\n\n")

if (!file.exists(annotated_vcf)) {
  stop("Annotated VCF not found (run 02_annotate_sQTL_snpeff.sh): ", annotated_vcf)
}

# --- Parse the annotated VCF --------------------------------------------------
cat("Reading annotated VCF file...\n")
annotated_variants <- parse_snpeff_vcf(annotated_vcf)
cat(sprintf("Extracted %d annotation entries for %d variants\n",
            nrow(annotated_variants),
            length(unique(annotated_variants$variant_id))))

# --- Join onto the significant sQTL ------------------------------------------
cat("\nReading sQTL results...\n")
sqtl <- fread(sqtl_file, header = TRUE)
sqtl_sig <- sqtl %>%
  filter(FDR < FDR_THRESHOLD) %>%
  # The VCF carries no "chr" prefix, so the join key is built the same way here
  mutate(variant_id = paste(gsub("^chr", "", CHR), POS, Allele1, Allele2,
                            sep = ":"))
cat(sprintf("Significant sQTLs: %d\n", nrow(sqtl_sig)))

cat("Merging annotations with sQTL data...\n")
sqtl_annotated <- sqtl_sig %>%
  left_join(annotated_variants, by = "variant_id")

annotated_output <- file.path(output_dir, "sQTL_with_snpeff_annotations.tsv")
fwrite(sqtl_annotated, annotated_output, sep = "\t", quote = FALSE, na = "NA")

cat(sprintf("\nAnnotation entries written: %d\n", nrow(sqtl_annotated)))
cat(sprintf("sQTL without any annotation: %d\n",
            sum(is.na(sqtl_annotated$Annotation))))
cat("\n=== DONE ===\n")
cat(sprintf("Annotated sQTLs saved to: %s\n", annotated_output))
