#!/usr/bin/env Rscript
# Enrichment per annotation type.
#
# Quantify how the functional annotations of the lead sQTL variants differ
# from those of the non-significant background variants (steps 04-05).
# SnpEff annotations are kept at full granularity, compound annotations
# included, and every annotation type is tested separately: for each type a
# 2x2 contingency table of presence/absence in lead versus background
# variants is tested with Fisher's exact test (two-sided). Effect sizes are
# odds ratios with 95% confidence intervals; p-values are adjusted across the
# tested annotation types with Benjamini-Hochberg.
#
# The log odds ratios are the sQTL weights shown in the forest plots:
# positive means enrichment in sQTL variants, negative means depletion.

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(ggplot2)
})

source("enrichment_functions.R")

# --- Paths (relative to this script's directory) -----------------------------
snpeff_dir <- "../../data/05_sQTL_mapping/SNPeff"
qtlbhm_dir <- "../../data/05_sQTL_mapping/qtlBHM"

sig_annotated_file <- file.path(snpeff_dir, "sQTL_with_snpeff_annotations.tsv")
bg_annotated_vcf   <- file.path(qtlbhm_dir, "background_variants_annotated.vcf")

if (!file.exists(sig_annotated_file)) {
  stop("Significant sQTL annotations not found: ", sig_annotated_file)
}
if (!file.exists(bg_annotated_vcf)) {
  stop("Background annotations not found (run step 05): ", bg_annotated_vcf)
}

# Annotation types seen fewer than this many times in both sets are skipped,
# because their odds ratios are unstable.
MIN_COUNT <- 10

cat("=== Enrichment analysis: all SnpEff annotation types ===\n\n")

# ============================================================
# LEAD sQTL VARIANTS
# ============================================================

cat("Reading significant sQTL annotations...\n")
sig_sqtls <- fread(sig_annotated_file, header = TRUE)
cat(sprintf("Significant sQTL entries: %d\n", nrow(sig_sqtls)))

sig_clean <- sig_sqtls %>%
  filter(!is.na(variant_id) & variant_id != "") %>%
  filter(!is.na(Annotation) & Annotation != "")
cat(sprintf("Entries with valid annotations: %d\n", nrow(sig_clean)))

# One lead variant per sQTL: the one with the smallest FDR
sig_lead_variants <- sig_clean %>%
  group_by(sQTL) %>%
  slice_min(FDR, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  pull(variant_id) %>%
  unique()
cat(sprintf("Unique lead sQTL variants: %d\n", length(sig_lead_variants)))

# ============================================================
# BACKGROUND VARIANTS
# ============================================================

cat("\nParsing background variant annotations from VCF...\n")
bg_annotated <- parse_snpeff_vcf(bg_annotated_vcf)
cat(sprintf("Background annotation entries: %d\n", nrow(bg_annotated)))

# ============================================================
# ANNOTATION TYPES
# ============================================================

cat("\n=== Identifying all annotation types ===\n")

all_annotation_types <- union(unique(sig_clean$Annotation),
                              unique(bg_annotated$Annotation))
cat(sprintf("Total unique annotation types: %d\n", length(all_annotation_types)))

cat("\nTop 20 annotation types in significant sQTLs:\n")
print(sig_clean %>% count(Annotation, name = "n_sig") %>%
        arrange(desc(n_sig)) %>% head(20))

cat("\nTop 20 annotation types in background:\n")
print(bg_annotated %>% count(Annotation, name = "n_bg") %>%
        arrange(desc(n_bg)) %>% head(20))

# Variant -> annotation assignments, exactly as SnpEff reported them
sig_variant_annotations <- sig_clean %>%
  filter(variant_id %in% sig_lead_variants) %>%
  select(variant_id, Annotation) %>%
  distinct()
bg_variant_annotations <- bg_annotated %>%
  select(variant_id, Annotation) %>%
  distinct()

sig_variant_ids <- unique(sig_variant_annotations$variant_id)
bg_variant_ids  <- unique(bg_variant_annotations$variant_id)

# ============================================================
# ENRICHMENT ANALYSIS
# ============================================================

cat("\n=== Running enrichment analysis for all annotation types ===\n\n")
cat(sprintf("Lead variants: %d\n", length(sig_variant_ids)))
cat(sprintf("Background variants: %d\n\n", length(bg_variant_ids)))

enrichment_results <- list()

for (annot in all_annotation_types) {
  n_sig_with <- sig_variant_annotations %>%
    filter(Annotation == annot) %>% pull(variant_id) %>% unique() %>% length()
  n_bg_with <- bg_variant_annotations %>%
    filter(Annotation == annot) %>% pull(variant_id) %>% unique() %>% length()

  if (n_sig_with < MIN_COUNT && n_bg_with < MIN_COUNT) next

  enrichment_results[[annot]] <- enrichment_fisher(
    annot, n_sig_with, length(sig_variant_ids),
    n_bg_with, length(bg_variant_ids))
}

enrichment_df <- do.call(rbind, enrichment_results)
rownames(enrichment_df) <- NULL
enrichment_df <- enrichment_df %>% rename(Annotation = label)
enrichment_df$padj <- p.adjust(enrichment_df$pvalue, method = "BH")
enrichment_df <- drop_infinite(enrichment_df)

cat("\n=== Enrichment results ===\n\n")
cat(sprintf("Tested %d annotation types (at least %d occurrences)\n\n",
            nrow(enrichment_df), MIN_COUNT))

print(enrichment_df %>%
        arrange(padj, desc(abs(log_odds))) %>%
        mutate(odds_ratio = round(odds_ratio, 3),
               across(c(log_odds, log_ci_lower, log_ci_upper), ~ round(.x, 3)),
               pvalue_fmt = formatC(pvalue, format = "e", digits = 2),
               padj_fmt   = formatC(padj, format = "e", digits = 2)) %>%
        select(Annotation, n_sig, n_bg, odds_ratio, log_odds,
               log_ci_lower, log_ci_upper, pvalue_fmt, padj_fmt))

enrichment_output <- file.path(qtlbhm_dir, "enrichment_all_annotations.tsv")
fwrite(enrichment_df, enrichment_output, sep = "\t", quote = FALSE,
       row.names = FALSE)
cat(sprintf("\nEnrichment results saved to: %s\n", enrichment_output))

# ============================================================
# FOREST PLOTS
# ============================================================

cat("\n=== Creating forest plots ===\n")

add_significance <- function(df) {
  df %>% mutate(Significance = case_when(padj < 0.01 ~ "FDR < 0.01",
                                         padj < 0.05 ~ "FDR < 0.05",
                                         TRUE ~ "Not significant"),
                Significance = factor(Significance,
                                      levels = c("FDR < 0.01", "FDR < 0.05",
                                                 "Not significant")))
}

sig_colors <- c("FDR < 0.01" = "#D32F2F", "FDR < 0.05" = "#FF9800",
                "Not significant" = "#757575")

forest_theme <- theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        plot.subtitle = element_text(size = 10, hjust = 0.5),
        axis.title = element_text(size = 12, face = "bold"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "top")

# Annotation types ordered by their sQTL weight
plot_data <- enrichment_df %>%
  filter(!is.na(log_odds)) %>%
  arrange(log_odds) %>%
  mutate(Annotation_display = factor(Annotation, levels = Annotation)) %>%
  add_significance()

pdf(file.path(qtlbhm_dir, "enrichment_forest_plot_all_annotations.pdf"),
    width = 10, height = max(6, nrow(plot_data) * 0.3))
print(
  ggplot(plot_data, aes(x = log_odds, y = Annotation_display,
                        color = Significance)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50",
               linewidth = 0.8) +
    geom_errorbar(aes(xmin = log_ci_lower, xmax = log_ci_upper),
                  width = 0.3, linewidth = 0.6) +
    geom_point(size = 3) +
    scale_color_manual(values = sig_colors) +
    labs(title = "sQTL Variant Enrichment Across SnpEff Annotation Types",
         subtitle = sprintf(
           "Lead sQTL variants vs background (n_sig=%d, n_bg=%d, min_count=%d)",
           length(sig_variant_ids), length(bg_variant_ids), MIN_COUNT),
         x = "sQTL Weight (Log Odds Ratio)", y = "SnpEff Annotation Type",
         color = "Significance") +
    forest_theme +
    theme(axis.text.y = element_text(size = 8),
          axis.text.x = element_text(size = 10)))
dev.off()

# Significant annotation types only
plot_data_sig <- plot_data %>% filter(padj < 0.05)

if (nrow(plot_data_sig) > 0) {
  pdf(file.path(qtlbhm_dir, "enrichment_forest_plot_significant_only.pdf"),
      width = 10, height = max(4, nrow(plot_data_sig) * 0.4))
  print(
    ggplot(plot_data_sig, aes(x = log_odds, y = Annotation_display,
                              color = Significance)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey50",
                 linewidth = 0.8) +
      geom_errorbar(aes(xmin = log_ci_lower, xmax = log_ci_upper),
                    width = 0.3, linewidth = 0.8) +
      geom_point(size = 4) +
      scale_color_manual(values = sig_colors) +
      labs(title = "Significantly Enriched/Depleted Annotation Types",
           subtitle = "FDR-adjusted p-value < 0.05",
           x = "sQTL Weight (Log Odds Ratio)", y = "SnpEff Annotation Type",
           color = "Significance") +
      forest_theme + theme(axis.text = element_text(size = 11)))
  dev.off()
  cat(sprintf("Significant annotation types (FDR < 0.05): %d\n",
              nrow(plot_data_sig)))
} else {
  cat("No significant enrichments found (FDR < 0.05)\n")
}

# The 20 strongest effects in either direction
plot_data_top <- enrichment_df %>%
  filter(!is.na(log_odds)) %>%
  arrange(desc(abs(log_odds))) %>%
  head(20) %>%
  arrange(log_odds) %>%
  mutate(Annotation_display = factor(Annotation, levels = Annotation)) %>%
  add_significance()

pdf(file.path(qtlbhm_dir, "enrichment_forest_plot_top20_effect_size.pdf"),
    width = 10, height = 8)
print(
  ggplot(plot_data_top, aes(x = log_odds, y = Annotation_display,
                            color = Significance)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50",
               linewidth = 0.8) +
    geom_errorbar(aes(xmin = log_ci_lower, xmax = log_ci_upper),
                  width = 0.3, linewidth = 0.8) +
    geom_point(size = 4) +
    scale_color_manual(values = sig_colors) +
    labs(title = "Top 20 Annotation Types by Effect Size",
         subtitle = "Ranked by absolute log odds ratio",
         x = "sQTL Weight (Log Odds Ratio)", y = "SnpEff Annotation Type",
         color = "Significance") +
    forest_theme + theme(axis.text = element_text(size = 11)))
dev.off()

cat("\n=== DONE ===\n")
cat(sprintf("  Annotation types tested: %d\n", nrow(enrichment_df)))
cat(sprintf("  Significant (FDR < 0.05): %d\n",
            sum(enrichment_df$padj < 0.05, na.rm = TRUE)))
cat(sprintf("Output written to %s\n", qtlbhm_dir))
