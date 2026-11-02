#!/usr/bin/env Rscript
# 03.GSEA_barplot_committed_aspc.R
# GSEA interaction term barplot for Committed_preadipocyte
# x = NES, y = pathway, color = p.adj, facet = GOBP / Hallmark

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(stringr)
})

# Project root is 7 levels up from this script
args0 <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args0[grep("--file=", args0)])
project_root <- normalizePath(file.path(dirname(normalizePath(script_path)),
                                        paste(rep("..", 7), collapse = "/")))

GSEA_DIR <- file.path(project_root,
                      "data/21.differential_expressed_genes/01.pseudobulk",
                      "04.only_VAT_vs_SAT_pseudobulk/linear_model/interation_term/tables/GSEA")
OUT_DIR  <- file.path(project_root,
                      "data/21.differential_expressed_genes/01.pseudobulk",
                      "04.only_VAT_vs_SAT_pseudobulk/linear_model/interation_term/plots")

CT <- "Committed_preadipocyte"

# ─── Load ─────────────────────────────────────────────────────────────────────
gobp <- as.data.frame(fread(file.path(GSEA_DIR,
  paste0(CT, ".GSEA.sig.GOBP_interaction.tsv"))))
hall <- as.data.frame(fread(file.path(GSEA_DIR,
  paste0(CT, ".GSEA.sig.Hallmark_interaction.tsv"))))

gobp$gene_set <- "GOBP"
hall$gene_set <- "Hallmark"

df <- bind_rows(gobp, hall)
cat(sprintf("GOBP: %d pathways, Hallmark: %d pathways\n", nrow(gobp), nrow(hall)))

# ─── Keep only 5 selected pathways ───────────────────────────────────────────
SELECTED <- c(
  "GOBP_EXTERNAL_ENCAPSULATING_STRUCTURE_ORGANIZATION",
  "HALLMARK_ADIPOGENESIS",
  "HALLMARK_MTORC1_SIGNALING"
)
df <- df %>% filter(pathway %in% SELECTED)
cat(sprintf("After filtering: %d pathways\n", nrow(df)))

# ─── Clean pathway names ──────────────────────────────────────────────────────
clean_pw_name <- function(x, wrap = 28) {
  x <- sub("^GOBP_", "", x)
  x <- sub("^HALLMARK_", "", x)
  x <- gsub("_", " ", x)
  x <- tolower(x)
  x <- sub("^(\\w)", "\\U\\1", x, perl = TRUE)
  stringr::str_wrap(x, width = wrap)
}
df$pathway_clean <- clean_pw_name(df$pathway)
df$pathway_clean <- gsub("organization", "org.", df$pathway_clean)

# order pathways by NES within each gene_set
df <- df %>%
  arrange(gene_set, NES) %>%
  mutate(
    pathway_clean  = factor(pathway_clean, levels = unique(pathway_clean)),
    neg_log10_padj = -log10(padj)
  )

DIV_COLORS <- c(
  "#543005","#8c510a","#bf812d","#dfc27d","#f6e8c3",
  "#f5f5f5",
  "#c7eae5","#80cdc1","#35978f","#01665e","#003c30"
)

make_barplot <- function(dat) {
  # scale limits computed from this dataset only
  nes_lim  <- max(abs(dat$NES)) * 1.05
  padj_lo  <- floor(min(dat$neg_log10_padj)  * 10) / 10
  padj_hi  <- ceiling(max(dat$neg_log10_padj) * 10) / 10

  ggplot(dat, aes(x = NES, y = pathway_clean, fill = neg_log10_padj)) +
  geom_col(width = 0.55) +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey40") +
  scale_fill_gradientn(
    colors = DIV_COLORS,
    limits = c(padj_lo, padj_hi),
    breaks = round(c(padj_lo, (padj_lo + padj_hi) / 2, padj_hi), 1),
    name   = expression(-log[10](adj.p)),
    guide  = guide_colorbar(
      title.position = "top",
      barwidth  = unit(2, "cm"),
      barheight = unit(0.22, "cm"),
      ticks     = FALSE
    )
  ) +
  scale_x_continuous(
    limits = c(-nes_lim, nes_lim),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  facet_wrap(~ gene_set, scales = "free_y", ncol = 1) +
  labs(x = "NES", y = NULL, title = paste0(CT, " - GSEA interaction")) +
  theme_classic(base_size = 6) +
  theme(
    plot.title       = element_text(size = 6, face = "bold", hjust = 0.5),
    strip.text       = element_text(size = 5.5, face = "bold"),
    strip.background = element_blank(),
    axis.text.y      = element_text(size = 5),
    axis.text.x      = element_text(size = 5),
    axis.title.x     = element_text(size = 5.5),
    legend.position  = "bottom",
    legend.direction = "horizontal",
    legend.title     = element_text(size = 5),
    legend.text      = element_text(size = 4.5),
    legend.margin    = margin(0, 0, 0, 0),
    legend.box.margin = margin(-4, 0, 0, 0),
    panel.spacing    = unit(0.3, "lines"),
    plot.margin      = margin(3, 3, 2, 3)
  )
}

# ─── top3 (main figure) ───────────────────────────────────────────────────────
p_top3 <- make_barplot(df)
out <- file.path(OUT_DIR, paste0(CT, "_GSEA_interaction_barplot_top3.pdf"))
ggsave(out, p_top3, width = 3.8, height = 1.8, device = "pdf")
cat(sprintf("Saved: %s (3.8 x 1.8 in)\n", basename(out)))

# ─── all significant pathways (supplementary) ────────────────────────────────
df_all <- bind_rows(
  gobp %>% mutate(gene_set = "GOBP"),
  hall %>% mutate(gene_set = "Hallmark")
) %>%
  mutate(pathway_clean = clean_pw_name(pathway, wrap = 35)) %>%
  arrange(gene_set, NES) %>%
  mutate(
    pathway_clean  = factor(pathway_clean, levels = unique(pathway_clean)),
    neg_log10_padj = -log10(padj)
  )

padj_lim <- c(0, ceiling(max(df_all$neg_log10_padj) * 1.1))
nes_lim  <- max(abs(df_all$NES)) * 1.05

p_all <- make_barplot(df_all)
n_paths <- nrow(df_all)
out_all <- file.path(OUT_DIR, paste0(CT, "_GSEA_interaction_barplot_all.pdf"))
ggsave(out_all, p_all, width = 4.5, height = 0.28 * n_paths + 1.2, device = "pdf")
cat(sprintf("Saved: %s (4.5 x %.1f in)\n", basename(out_all),
            0.28 * n_paths + 1.2))
