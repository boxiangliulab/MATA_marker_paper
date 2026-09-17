#!/usr/bin/env Rscript
# Significant sQTL variants as VCF.
#
# Write the unique variants of the significant sQTL (FDR < 0.05) as a VCF so
# that SnpEff can annotate them in step 02. Chromosome names are written
# without the "chr" prefix, which is what the GRCh38 SnpEff database expects,
# and the records are sorted by chromosome and position.

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
})

source("enrichment_functions.R")

# --- Paths (relative to this script's directory) -----------------------------
results_dir <- "../../data/05_sQTL_mapping/results"
snpeff_dir  <- "../../data/05_sQTL_mapping/SNPeff"

dir.create(snpeff_dir, recursive = TRUE, showWarnings = FALSE)

sqtl_file <- file.path(results_dir, "All_sQTL_SAT_and_VAT_with_beta.tsv")
vcf_file <- file.path(snpeff_dir, "sQTL_variants_FDR05.vcf")

FDR_THRESHOLD <- 0.05

cat("=== Convert sQTL results to VCF format ===\n\n")

if (!file.exists(sqtl_file)) {
  stop("sQTL summary not found (run 01_mapping/03_summarize_sQTL.R): ", sqtl_file)
}

cat("Reading sQTL results...\n")
sqtl <- fread(sqtl_file, header = TRUE)
cat(sprintf("Total sQTLs: %d\n", nrow(sqtl)))

cat(sprintf("Filtering for FDR < %.2f...\n", FDR_THRESHOLD))
sqtl_sig <- sqtl %>% filter(FDR < FDR_THRESHOLD)
cat(sprintf("Significant sQTLs: %d\n", nrow(sqtl_sig)))

# One record per variant, not per association
cat("Extracting unique variants...\n")
variants <- sqtl_sig %>%
  select(CHR, POS, MarkerID, Allele1, Allele2) %>%
  distinct() %>%
  arrange(CHR, POS)
cat(sprintf("Unique significant variants: %d\n", nrow(variants)))

variants <- variants %>%
  mutate(CHROM  = gsub("^chr", "", CHR),   # SnpEff GRCh38 has no "chr" prefix
         ID     = ".",                     # no rsID available here
         REF    = Allele1,
         ALT    = Allele2,
         QUAL   = ".",
         FILTER = "PASS",
         INFO   = ".") %>%
  select(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO)

# Sort numerically, with the sex and mitochondrial chromosomes last
variants <- variants %>%
  mutate(CHR_order = case_when(
    CHROM %in% as.character(1:22) ~ as.integer(CHROM),
    CHROM == "X"  ~ 23L,
    CHROM == "Y"  ~ 24L,
    CHROM == "M"  ~ 25L,
    CHROM == "MT" ~ 25L,
    TRUE ~ 99L)) %>%
  arrange(CHR_order, POS) %>%
  select(-CHR_order)

cat(sprintf("\nWriting VCF file to: %s\n", vcf_file))

vcf_header <- c(
  "##fileformat=VCFv4.2",
  sprintf("##fileDate=%s", format(Sys.Date(), "%Y%m%d")),
  "##source=sQTL_analysis",
  "##reference=GRCh38",
  paste0("##contig=<ID=", unique(variants$CHROM), ">"),
  "##INFO=<ID=.,Number=0,Type=Flag,Description=\"No additional info\">",
  "##FILTER=<ID=PASS,Description=\"All filters passed\">",
  paste0("#", paste(c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
                      "INFO"), collapse = "\t"))
)

writeLines(vcf_header, vcf_file)
fwrite(variants, vcf_file, sep = "\t", quote = FALSE,
       row.names = FALSE, col.names = FALSE, append = TRUE)

cat(sprintf("VCF file created with %d variants\n", nrow(variants)))

cat("\nVariants per chromosome:\n")
print(variants %>%
        group_by(CHROM) %>%
        summarise(n_variants = n(), .groups = "drop") %>%
        arrange(as.integer(ifelse(CHROM %in% c("X", "Y", "M", "MT"), "99",
                                  CHROM))))

cat("\n=== DONE ===\n")
cat(sprintf("VCF file saved to: %s\n", vcf_file))
