#!/usr/bin/env Rscript
# Background variant set.
#
# Build the background against which the functional enrichment of sQTL
# variants is measured: variants that were tested for association but are not
# significant. Up to 1,000 variants with p > 0.01 are sampled from every
# cell-type/depot result file, variants that are significant anywhere are
# removed, and the remainder is written as a VCF ready for SnpEff (step 05).

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

source("enrichment_functions.R")

# --- Paths (relative to this script's directory) -----------------------------
results_dir <- "../../data/05_sQTL_mapping/results"
qtlbhm_dir  <- "../../data/05_sQTL_mapping/qtlBHM"

dir.create(qtlbhm_dir, recursive = TRUE, showWarnings = FALSE)

# Number of non-significant variants sampled per result file
N_SAMPLE_PER_FILE <- 1000
# Variants with a nominal p-value above this are considered non-significant
P_THRESHOLD <- 0.01

cat("=== Creating Background Variant Set ===\n\n")

# --- Significant variants (to be excluded from the background) ---------------
cat("Reading significant sQTL data...\n")
sig_sqtls <- fread(file.path(results_dir,
                             "All_sQTL_SAT_and_VAT_with_beta.tsv"))
cat(sprintf("Significant sQTLs: %d\n", nrow(sig_sqtls)))

sig_variants <- sig_sqtls %>%
  mutate(variant_id = paste(gsub("^chr", "", CHR), POS, Allele1, Allele2,
                            sep = ":")) %>%
  select(variant_id, CHR, POS, Allele1, Allele2) %>%
  distinct()

cat(sprintf("Unique significant variants: %d\n", nrow(sig_variants)))

# --- Sample non-significant variants from every result file ------------------
cat("\nSampling background variants from full QTL results...\n")

result_files <- list.files(results_dir, pattern = "\\.result$",
                           full.names = TRUE)
cat(sprintf("Found %d full result files\n", length(result_files)))

background_variants_list <- list()

for (file_path in result_files) {
  cat(sprintf("  Processing: %s\n", basename(file_path)))

  # ISSAC output: splice_site, snp, pvalue, beta, se (no header)
  qtl_data <- fread(file_path, header = FALSE,
                    col.names = c("splice_site", "snp", "pvalue", "beta", "se"))
  cat(sprintf("    Total tests: %d\n", nrow(qtl_data)))

  background_subset <- qtl_data %>%
    filter(pvalue > P_THRESHOLD) %>%
    sample_n(min(N_SAMPLE_PER_FILE, n())) %>%
    select(snp) %>%
    distinct()

  cat(sprintf("    Sampled background variants: %d\n", nrow(background_subset)))
  background_variants_list[[length(background_variants_list) + 1]] <- background_subset
}

background_variants_raw <- do.call(rbind, background_variants_list) %>% distinct()
cat(sprintf("\nTotal background variant IDs: %d\n", nrow(background_variants_raw)))

# --- Parse the coordinates out of the variant IDs ----------------------------
cat("\nParsing SNP coordinates...\n")
background_variants <- background_variants_raw %>%
  mutate(
    CHR     = sub("^(chr[^:]+):.*", "\\1", snp),
    POS     = as.integer(sub("^chr[^:]+:(\\d+):.*", "\\1", snp)),
    Allele1 = sub("^chr[^:]+:\\d+:([^:]+):.*", "\\1", snp),
    Allele2 = sub("^chr[^:]+:\\d+:[^:]+:(.+)$", "\\1", snp),
    variant_id = paste(gsub("^chr", "", CHR), POS, Allele1, Allele2, sep = ":")
  ) %>%
  select(variant_id, CHR, POS, Allele1, Allele2) %>%
  distinct()

background_variants_clean <- background_variants %>%
  anti_join(sig_variants, by = "variant_id")

cat(sprintf("Background variants (excluding significant): %d\n",
            nrow(background_variants_clean)))

# --- Write the VCF for SnpEff -------------------------------------------------
cat("\nCreating VCF file for background variants...\n")
vcf_file <- file.path(qtlbhm_dir, "background_variants.vcf")

vcf_data <- background_variants_clean %>%
  mutate(CHROM = gsub("^chr", "", CHR), ID = ".", REF = Allele1, ALT = Allele2,
         QUAL = ".", FILTER = "PASS", INFO = ".") %>%
  mutate(CHR_order = case_when(
    CHROM %in% as.character(1:22) ~ as.integer(CHROM),
    CHROM == "X"  ~ 23L,
    CHROM == "Y"  ~ 24L,
    CHROM == "M"  ~ 25L,
    CHROM == "MT" ~ 25L,
    TRUE ~ 99L)) %>%
  arrange(CHR_order, POS) %>%
  select(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO)

vcf_header <- c(
  "##fileformat=VCFv4.2",
  sprintf("##fileDate=%s", format(Sys.Date(), "%Y%m%d")),
  "##source=background_variants",
  "##reference=GRCh38",
  paste0("##contig=<ID=", unique(vcf_data$CHROM), ">"),
  "##INFO=<ID=.,Number=0,Type=Flag,Description=\"No additional info\">",
  "##FILTER=<ID=PASS,Description=\"All filters passed\">",
  paste0("#", paste(c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
                      "INFO"), collapse = "\t"))
)

writeLines(vcf_header, vcf_file)
fwrite(vcf_data, vcf_file, sep = "\t", quote = FALSE,
       row.names = FALSE, col.names = FALSE, append = TRUE)

cat(sprintf("Background VCF created: %s\n", vcf_file))
cat(sprintf("Background variants: %d\n", nrow(vcf_data)))

background_list_file <- file.path(qtlbhm_dir, "background_variant_ids.txt")
writeLines(background_variants_clean$variant_id, background_list_file)
cat(sprintf("Background variant IDs saved: %s\n", background_list_file))

cat("\n=== DONE ===\n")
cat("Next step: annotate the VCF with SnpEff (05_annotate_background_snpeff.sh)\n")
