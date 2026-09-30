library(ggplot2)
library(patchwork)

# Project root is 8 levels up from this script
args0        <- commandArgs(trailingOnly = FALSE)
script_path  <- sub("--file=", "", args0[grep("--file=", args0)])
project_root <- normalizePath(file.path(dirname(normalizePath(script_path)),
                                        paste(rep("..", 8), collapse = "/")))

tss_col <- "#80b1d3"
tes_col <- "#fb8072"
df <- read.table(file.path(project_root, "data/09.eQTL/03.Saige_QTL/02.output/05.tss_region/ALL_leadSNP_distance_to_TSS_and_TES_all_celltype.tsv"), header=T, sep="\t")
# data.frame with columns: distance_tss, distance_tes (in bp)
df$distance_tss_kb <- df$distance_tss / 1000
df$distance_tes_kb <- df$distance_tes / 1000

# reshape to long format for ggplot
df_long <- rbind(
  data.frame(type = "TSS", dist_kb = df$distance_tss_kb),
  data.frame(type = "TES", dist_kb = df$distance_tes_kb)
)
df_long <- df_long[is.finite(df_long$dist_kb), ]

cols <- c("TSS" = "#80b1d3", "TES" = "#fb8072")

# ---------------------------
# 1) Separate plots (stacked)
# ---------------------------
p_tss <- ggplot(data.frame(x = df$distance_tss_kb), aes(x = x)) +
  geom_density(fill = tss_col, color = "black", linewidth = 0.4, alpha = 0.55) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", linewidth = 0.8) +
  labs(x = "Relative Distance to TSS (kb)", y = "Density") +
  theme_classic(base_size = 14) +
  theme(axis.title = element_text(face = "bold")) + coord_cartesian(xlim = c(-500, 500))

p_tes <- ggplot(data.frame(x = df$distance_tes_kb), aes(x = x)) +
  geom_density(fill = tes_col, color = "black", linewidth = 0.4, alpha = 0.55) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", linewidth = 0.8) +
  labs(x = "Relative Distance to TES (kb)", y = "Density") +
  theme_classic(base_size = 14) +
  theme(axis.title = element_text(face = "bold")) + coord_cartesian(xlim = c(-500, 500))

p_sep <- p_tss / p_tes
ggsave(file.path(project_root, "data/09.eQTL/03.Saige_QTL/02.output/05.tss_region/TSS_TES_density_separate.pdf"), p_sep, width = 6.2, height = 9)

# ---------------------------
# 2) Overlay plot (together)
# ---------------------------
df_long <- rbind(
  data.frame(type = "TSS", dist_kb = df$distance_tss_kb),
  data.frame(type = "TES", dist_kb = df$distance_tes_kb)
)
df_long <- df_long[is.finite(df_long$dist_kb), ]

fill_cols <- c("TSS" = tss_col, "TES" = tes_col)

p_overlay <- ggplot(df_long, aes(x = dist_kb, fill = type)) +
  geom_density(aes(group = type), color = "black", linewidth = 0.4, alpha = 0.35, adjust = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", linewidth = 0.8) +
  scale_fill_manual(values = fill_cols) +
  labs(x = "Relative distance (kb)", y = "Density") +
  theme_classic(base_size = 14) +
  theme(
    axis.title = element_text(face = "bold"),
    legend.title = element_blank()
  ) + coord_cartesian(xlim = c(-500, 500))

ggsave(file.path(project_root, "data/09.eQTL/03.Saige_QTL/02.output/05.tss_region/TSS_TES_density_overlay.pdf"), p_overlay, width = 6.5, height = 4.5)

