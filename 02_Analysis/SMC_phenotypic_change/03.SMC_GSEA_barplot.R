library(data.table)
library(dplyr)
library(ggplot2)
library(stringr)

# ============================================================
# Paths
# ============================================================
# Project root is 7 levels up from this script
args0 <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args0[grep("--file=", args0)])
project_root <- normalizePath(file.path(dirname(normalizePath(script_path)),
                                        paste(rep("..", 7), collapse = "/")))

gsea_file <- file.path(project_root, "data/21.differential_expressed_genes/01.pseudobulk/04.only_VAT_vs_SAT_pseudobulk/linear_model/interation_term/tables/GSEA/SMC.GSEA.sig.GOBP_interaction.tsv")
out_dir   <- file.path(project_root, "data/35.vascular/01.DEG")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# Read & clean
# ============================================================
df <- fread(gsea_file, data.table = FALSE)

DIV_COLORS <- c(
  "#543005","#8c510a","#bf812d","#dfc27d","#f6e8c3",
  "#f5f5f5",
  "#c7eae5","#80cdc1","#35978f","#01665e","#003c30"
)

df <- df %>%
  arrange(desc(abs(NES))) %>%
  slice_head(n = 30) %>%
  mutate(
    pathway_clean = pathway %>%
      str_remove("^GOBP_") %>%
      str_replace_all("_", " ") %>%
      str_to_title() %>%
      str_wrap(width = 45),
    neglog10_padj = -log10(padj)
  ) %>%
  arrange(NES) %>%
  mutate(pathway_clean = factor(pathway_clean, levels = pathway_clean))

# ============================================================
# Plot
# ============================================================
x_lim     <- max(abs(df$NES), na.rm = TRUE) + 0.2
padj_lim  <- c(0, ceiling(max(df$neglog10_padj) * 1.1))

p <- ggplot(df, aes(x = NES, y = pathway_clean, fill = neglog10_padj)) +
  geom_bar(stat = "identity", width = 0.75) +
  geom_vline(xintercept = 0, linewidth = 0.4, color = "grey30") +
  scale_fill_gradientn(
    colors = DIV_COLORS,
    limits = padj_lim,
    name   = expression(-log[10](adj.p)),
    guide  = guide_colorbar(
      title.position = "top",
      barwidth  = unit(3, "cm"),
      barheight = unit(0.35, "cm"),
      ticks     = FALSE
    )
  ) +
  coord_cartesian(xlim = c(-x_lim, x_lim), clip = "off") +
  labs(
    title = "SMC: Depot x BMI interaction - GO Biological Process (FDR < 0.05)",
    x = "Normalized Enrichment Score (NES)",
    y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.y        = element_text(size = 8),
    legend.position    = "top",
    plot.title         = element_text(size = 11, face = "bold")
  )

# ============================================================
# Save
# ============================================================
out_file <- file.path(out_dir, "SMC.GSEA.GOBP_interaction.top30.barplot.pdf")
if (file.exists(out_file)) file.remove(out_file)
ggsave(out_file, p, width = 4.5, height = 6, device = "pdf")

cat("Saved:", out_file, "\n")
