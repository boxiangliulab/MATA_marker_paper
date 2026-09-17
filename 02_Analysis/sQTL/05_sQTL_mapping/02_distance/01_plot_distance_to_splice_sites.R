#!/usr/bin/env Rscript
# Distance from sQTL to splice sites.
#
# Describe where significant sQTL variants sit relative to the splice site
# they affect: summary statistics, the fraction of variants within 1 kb /
# 10 kb / 100 kb / 1 Mb, and density plots of the signed and absolute
# distance (overall, zoomed to +/-100 kb, and split by depot).

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(scales)
})

cat("=== Distance from SNPs to Splice Sites Analysis ===\n\n")

# --- Paths (relative to this script's directory) -----------------------------
results_dir <- "../../data/05_sQTL_mapping/results"
output_dir  <- "../../data/05_sQTL_mapping/distance"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# --- Read the combined sQTL table -------------------------------------------
cat("Reading sQTL data...\n")
sqtl_file <- file.path(results_dir, "All_sQTL_SAT_and_VAT_with_beta.tsv")
sqtls <- fread(sqtl_file, header = TRUE)

cat(sprintf("Total sQTL associations: %d\n", nrow(sqtls)))
cat(sprintf("Unique sQTLs: %d\n", length(unique(sqtls$sQTL))))

# --- Parse the splice-site coordinate ----------------------------------------
# Splice-site IDs look like chr15:-:34582749 (chromosome, strand, position).
cat("\nParsing splice site positions...\n")
sqtls$splice_chr    <- sub("^(chr[^:]+):.*", "\\1", sqtls$splice_site)
sqtls$splice_strand <- sub("^chr[^:]+:([^:]+):.*", "\\1", sqtls$splice_site)
sqtls$splice_pos    <- as.integer(sub("^chr[^:]+:[^:]+:(\\d+).*", "\\1",
                                      sqtls$splice_site))

# Signed distance: negative = variant upstream of the splice site
sqtls$distance <- sqtls$POS - sqtls$splice_pos

cat(sprintf("Successfully parsed: %d/%d\n",
            sum(!is.na(sqtls$distance)), nrow(sqtls)))

# --- Summary statistics ------------------------------------------------------
cat("\n=== Distance Statistics ===\n")
cat(sprintf("Mean distance: %.1f bp\n", mean(abs(sqtls$distance), na.rm = TRUE)))
cat(sprintf("Median distance: %.1f bp\n", median(abs(sqtls$distance), na.rm = TRUE)))
cat(sprintf("Min distance: %d bp\n", min(abs(sqtls$distance), na.rm = TRUE)))
cat(sprintf("Max distance: %d bp\n", max(abs(sqtls$distance), na.rm = TRUE)))

cat("\nDistance distribution:\n")
for (cutoff in c(1000, 10000, 100000, 1000000)) {
  cat(sprintf("  Within %s: %d (%.1f%%)\n",
              format(cutoff, big.mark = ","),
              sum(abs(sqtls$distance) <= cutoff, na.rm = TRUE),
              100 * mean(abs(sqtls$distance) <= cutoff, na.rm = TRUE)))
}

# --- Save the per-association distances --------------------------------------
distance_data <- sqtls %>%
  select(sQTL, CHR, POS, splice_site, splice_pos, distance, depot, celltype,
         gene, FDR) %>%
  arrange(abs(distance))

fwrite(distance_data, file.path(output_dir, "sQTL_distances_to_splice_sites.tsv"),
       sep = "\t", quote = FALSE)
cat(sprintf("\nDistance data saved: %s\n",
            file.path(output_dir, "sQTL_distances_to_splice_sites.tsv")))

# ============================================================
# DENSITY PLOTS
# ============================================================

cat("\n=== Creating density plots ===\n")

plot_data <- sqtls %>%
  filter(!is.na(distance)) %>%
  mutate(abs_distance = abs(distance),
         distance_kb = distance / 1000,
         abs_distance_kb = abs(distance) / 1000)

# Shared x-axis for the log-scale plots
log_breaks <- c(0.1, 1, 10, 100, 1000)
log_labels <- c("100 bp", "1 kb", "10 kb", "100 kb", "1 Mb")

# --- 1. Absolute distance, log scale -----------------------------------------
pdf(file.path(output_dir, "distance_density_overall.pdf"), width = 10, height = 6)
p1 <- ggplot(plot_data, aes(x = abs_distance_kb)) +
  geom_density(fill = "#2196F3", alpha = 0.6, color = "#1976D2", linewidth = 1) +
  scale_x_log10(breaks = log_breaks, labels = log_labels) +
  labs(title = "Distance from sQTL Variants to Splice Sites",
       x = "Distance (absolute value, log scale)", y = "Density") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 11),
        panel.grid.minor = element_blank())
print(p1)
dev.off()
cat("Saved: distance_density_overall.pdf\n")

# --- 2. Signed distance (direction relative to the splice site) --------------
pdf(file.path(output_dir, "distance_density_signed.pdf"), width = 6, height = 5)
p2 <- ggplot(plot_data, aes(x = distance_kb)) +
  geom_density(fill = "grey60", alpha = 0.6, color = "grey30", linewidth = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black",
             linewidth = 0.8) +
  scale_x_continuous(limits = c(-1000, 1000), breaks = seq(-1000, 1000, 250),
                     labels = function(x) paste0(x, " kb")) +
  labs(title = "Signed Distance from sQTL Variants to Splice Sites",
       subtitle = "Negative = upstream of splice site, Positive = downstream",
       x = "Distance (signed, kb)", y = "Density") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        plot.subtitle = element_text(size = 10, hjust = 0.5),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 11))
print(p2)
dev.off()
cat("Saved: distance_density_signed.pdf\n")

# --- 3. Zoomed to +/-100 kb ---------------------------------------------------
pdf(file.path(output_dir, "distance_density_zoomed.pdf"), width = 10, height = 6)
plot_data_zoom <- plot_data %>% filter(abs_distance_kb <= 100)
p3 <- ggplot(plot_data_zoom, aes(x = distance_kb)) +
  geom_density(fill = "#FF9800", alpha = 0.6, color = "#F57C00", linewidth = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red",
             linewidth = 0.8) +
  scale_x_continuous(limits = c(-100, 100), breaks = seq(-100, 100, 20),
                     labels = function(x) paste0(x, " kb")) +
  labs(title = "Distance from sQTL Variants to Splice Sites (<=100 kb)",
       subtitle = sprintf("n = %s associations (%.1f%% of total)",
                          format(nrow(plot_data_zoom), big.mark = ","),
                          100 * nrow(plot_data_zoom) / nrow(plot_data)),
       x = "Distance (signed, kb)", y = "Density") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        plot.subtitle = element_text(size = 10, hjust = 0.5),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 11),
        panel.grid.minor = element_blank())
print(p3)
dev.off()
cat("Saved: distance_density_zoomed.pdf\n")

# --- 4. By depot --------------------------------------------------------------
pdf(file.path(output_dir, "distance_density_by_depot.pdf"), width = 10, height = 6)
p4 <- ggplot(plot_data, aes(x = abs_distance_kb, fill = depot)) +
  geom_density(alpha = 0.5, linewidth = 0.8) +
  scale_fill_manual(values = c("SAT" = "#E91E63", "VAT" = "#9C27B0")) +
  scale_x_log10(breaks = log_breaks, labels = log_labels) +
  labs(title = "Distance from sQTL Variants to Splice Sites by Depot",
       x = "Distance (absolute value, log scale)", y = "Density",
       fill = "Depot") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 11),
        panel.grid.minor = element_blank(),
        legend.position = "top")
print(p4)
dev.off()
cat("Saved: distance_density_by_depot.pdf\n")

# --- 5. Histogram with log-spaced bins ---------------------------------------
pdf(file.path(output_dir, "distance_histogram_log_bins.pdf"), width = 10, height = 6)
p5 <- ggplot(plot_data, aes(x = abs_distance_kb)) +
  geom_histogram(bins = 50, fill = "#673AB7", color = "white", alpha = 0.8) +
  scale_x_log10(breaks = log_breaks, labels = log_labels) +
  labs(title = "Distribution of sQTL Variant Distances to Splice Sites",
       x = "Distance (absolute value, log scale)", y = "Count") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 11),
        panel.grid.minor = element_blank())
print(p5)
dev.off()
cat("Saved: distance_histogram_log_bins.pdf\n")

cat("\n=== DONE ===\n")
cat(sprintf("All output written to %s\n", output_dir))
