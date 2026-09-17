#!/usr/bin/env Rscript
# Two-level multiple-testing correction.
#

suppressPackageStartupMessages(library(dplyr))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript multiple_testing_correction.R <merged_result> <out_lead> <out_sig>")
}

merged_file <- args[1]
out_lead    <- args[2]
out_sig     <- args[3]

FDR_THRESHOLD <- 0.05

# The concatenated ISSAC output has no header.
df <- read.table(merged_file, header = FALSE)
colnames(df) <- c("splice_site", "SNP", "pvalue", "beta", "se")

# --- Level 1: Bonferroni within each splice site -----------------------------
df <- df %>%
  group_by(splice_site) %>%
  mutate(n_variants = n(),
         within_site_bonf = pmin(pvalue * n_variants, 1)) %>%
  ungroup()

# Lead variant per splice site (smallest corrected p-value)
df_min <- df %>%
  group_by(splice_site) %>%
  arrange(within_site_bonf) %>%
  slice(1) %>%
  ungroup()

# --- Level 2: Benjamini-Hochberg across splice sites -------------------------
df_min <- df_min %>%
  mutate(global_FDR = p.adjust(within_site_bonf, method = "BH")) %>%
  arrange(global_FDR, within_site_bonf)

df_min_sig <- df_min %>% filter(global_FDR < FDR_THRESHOLD)

write.table(df_min, file = out_lead, sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(df_min_sig, file = out_sig, sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = TRUE)

cat(sprintf("Total splice sites tested: %d\n", nrow(df_min)))
cat(sprintf("Significant sQTL (global FDR < %.2f): %d\n",
            FDR_THRESHOLD, nrow(df_min_sig)))
