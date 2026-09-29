#!/usr/bin/env Rscript
# 04.GSEA_enrichment_committed_aspc.R
# GSEA enrichment (running score + barcode) plots for 5 selected pathways
# Committed_preadipocyte interaction term — supplementary figure
#
# Ranked list: delta_logFC = logFC_obese - logFC_non_obese (VAT vs SAT edgeR)
#   positive = VAT-enriched effect with BMI
#   negative = SAT-enriched effect with BMI
# NES / padj annotated from original limma-based GSEA results (.sig tsv files)
#
# Barcode coloring:
#   orange (#e6550d) = fibrosis-related genes
#   blue   (#2171b5) = adipogenesis-related genes
#   grey             = other pathway genes

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(msigdbr)
  library(ggrepel)
})

# ─── Paths ────────────────────────────────────────────────────────────────────
BASE_DIR <- paste0("/data/projects/11003054/e1101943/project/01.adipose_sqtl/",
                   "data/21.differential_expressed_genes/01.pseudobulk/",
                   "04.only_VAT_vs_SAT_pseudobulk")
NON_ALL  <- file.path(BASE_DIR, "BMIgrp_non_obese/01.DEG/01.table",
                      "Committed_preadipocyte.non_obese.paired.edgeR.all.tsv")
OB_ALL   <- file.path(BASE_DIR, "BMIgrp_obese/01.DEG/01.table",
                      "Committed_preadipocyte.obese.paired.edgeR.all.tsv")
GSEA_DIR <- file.path(BASE_DIR,
                      "linear_model/interation_term/tables/GSEA")
OUT_DIR  <- file.path(BASE_DIR,
                      "linear_model/interation_term/plots")
CT <- "Committed_preadipocyte"

# ─── Selected pathways ────────────────────────────────────────────────────────
SELECTED <- c(
  "GOBP_EXTERNAL_ENCAPSULATING_STRUCTURE_ORGANIZATION",
  "HALLMARK_ADIPOGENESIS",
  "HALLMARK_MTORC1_SIGNALING"
)

# ─── Genes to highlight on barcodes ──────────────────────────────────────────
# Keep to ≤5 per category to avoid crowding
FIBROSIS_GENES <- c("COL1A1", "COL1A2", "POSTN", "LOX", "MMP14")
ADIPO_GENES    <- c("ADIPOR2", "MGLL", "ANGPTL4", "UCP2", "ESRRA")


HI_COLORS <- c(Fibrosis = "#e6550d", Adipogenesis = "#2171b5", Other = "grey78")
HI_LW     <- c(Fibrosis = 0.55, Adipogenesis = 0.55, Other = 0.18)

# ─── Build ranked list (delta logFC) ─────────────────────────────────────────
cat("Building ranked gene list (delta logFC) ...\n")
non <- fread(NON_ALL, select = c("gene", "logFC")) %>%
  rename(logFC_non = logFC)
ob  <- fread(OB_ALL,  select = c("gene", "logFC")) %>%
  rename(logFC_ob = logFC)

merged <- inner_join(as.data.frame(non), as.data.frame(ob), by = "gene") %>%
  mutate(delta = logFC_ob - logFC_non) %>%
  filter(is.finite(delta))

stat <- setNames(merged$delta, merged$gene)
stat <- sort(stat, decreasing = TRUE)
N    <- length(stat)
cat(sprintf("  Ranked genes: %d\n", N))

# ─── Load NES / padj from original GSEA results ───────────────────────────────
gobp_sig <- as.data.frame(fread(file.path(GSEA_DIR,
  paste0(CT, ".GSEA.sig.GOBP_interaction.tsv"))))
hall_sig <- as.data.frame(fread(file.path(GSEA_DIR,
  paste0(CT, ".GSEA.sig.Hallmark_interaction.tsv"))))
gsea_ann <- bind_rows(gobp_sig, hall_sig) %>%
  filter(pathway %in% SELECTED) %>%
  select(pathway, NES, padj)
cat("NES / padj loaded:\n")
print(gsea_ann)

# ─── Get pathway gene sets from msigdbr ──────────────────────────────────────
cat("Loading gene sets from msigdbr ...\n")
bp_sets <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")
h_sets  <- msigdbr(species = "Homo sapiens", category = "H")
pw_list <- split(
  c(bp_sets$gene_symbol, h_sets$gene_symbol),
  c(bp_sets$gs_name,     h_sets$gs_name)
)[SELECTED]

# ─── Running ES ───────────────────────────────────────────────────────────────
calc_es <- function(stats, genes) {
  n   <- length(stats)
  hit <- names(stats) %in% genes
  nh  <- sum(hit)
  if (nh == 0) return(data.frame(rank = seq_len(n), es = 0))
  wt  <- sum(abs(stats[hit]))
  step <- ifelse(hit, abs(stats) / wt, -1 / (n - nh))
  data.frame(rank = seq_len(n), es = cumsum(step))
}

# ─── Helper: clean pathway name ───────────────────────────────────────────────
clean_pw <- function(x) {
  x <- sub("^GOBP_", "", x)
  x <- sub("^HALLMARK_", "", x)
  x <- gsub("_", " ", x)
  x <- tolower(x)
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

# ─── Plot function ────────────────────────────────────────────────────────────
make_panel <- function(pw_name, blue_genes = ADIPO_GENES) {
  genes  <- pw_list[[pw_name]]
  ann    <- gsea_ann[gsea_ann$pathway == pw_name, ]
  nes    <- ann$NES
  padj   <- ann$padj

  # Running ES
  es_df <- calc_es(stat, genes)

  # Barcode positions
  idx <- which(names(stat) %in% genes)
  bc_df <- data.frame(
    rank = idx,
    gene = names(stat)[idx]
  ) %>%
    mutate(category = case_when(
      gene %in% FIBROSIS_GENES ~ "Fibrosis",
      gene %in% blue_genes      ~ "Adipogenesis",
      TRUE                     ~ "Other"
    ),
    category = factor(category, levels = c("Fibrosis", "Adipogenesis", "Other")))

  hi_df <- filter(bc_df, category != "Other")

  # Color / limits
  es_col  <- if (nes > 0) "#01665e" else "#8c510a"
  es_fill <- if (nes > 0) "#c7eae5" else "#dfc27d"
  y_peak  <- if (nes > 0) max(es_df$es) else min(es_df$es)
  y_lo    <- min(es_df$es, 0) * 1.2
  y_hi    <- max(es_df$es, 0) * 1.2
  if (y_hi == 0) y_hi <-  0.05
  if (y_lo == 0) y_lo <- -0.05

  padj_fmt <- formatC(padj, format = "e", digits = 1)

  # ── ES panel ──
  p_es <- ggplot(es_df, aes(x = rank, y = es)) +
    geom_area(fill = es_fill, alpha = 0.55) +
    geom_line(color = es_col, linewidth = 0.65) +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "grey45") +
    annotate("text",
      x = N * 0.97,
      y = y_hi * 0.91,
      label = paste0("NES = ", sprintf("%.2f", nes),
                     "\nadj.p = ", padj_fmt),
      hjust = 1, vjust = 1, size = 2.1, lineheight = 1.3, color = "grey15"
    ) +
    scale_x_continuous(expand = c(0, 0), limits = c(1, N)) +
    scale_y_continuous(
      limits = c(y_lo, y_hi),
      breaks = pretty(c(y_lo, y_hi), n = 3)
    ) +
    labs(title = clean_pw(pw_name), x = NULL, y = "Enrichment score") +
    theme_classic(base_size = 7.5) +
    theme(
      plot.title      = element_text(size = 7, face = "bold", hjust = 0.5,
                                     margin = margin(b = 1)),
      axis.text.x     = element_blank(),
      axis.ticks.x    = element_blank(),
      axis.line.x     = element_blank(),
      axis.text.y     = element_text(size = 5.5),
      axis.title.y    = element_text(size = 6),
      plot.margin     = margin(4, 3, 0, 3)
    )

  # ── Barcode panel ──
  p_bc <- ggplot(bc_df) +
    geom_segment(aes(x = rank, xend = rank, y = 0, yend = 1,
                     color = category, linewidth = category)) +
    scale_color_manual(values = HI_COLORS, guide = "none") +
    scale_linewidth_manual(values = HI_LW, guide = "none") +
    {if (nrow(hi_df) > 0)
      geom_text_repel(
        data    = hi_df,
        mapping = aes(x = rank, y = 1.02, label = gene, color = category),
        inherit.aes  = FALSE,
        size         = 1.8,
        angle        = 90,
        hjust        = 0,
        vjust        = 0.4,
        direction    = "x",
        nudge_y      = 0.08,
        segment.size = 0.2,
        segment.color = "grey60",
        max.overlaps  = 30,
        min.segment.length = 0
      )} +
    scale_x_continuous(
      expand = c(0, 0),
      limits = c(1, N),
      breaks = c(1, round(N / 2), N),
      labels = c("1", paste0(round(N / 2 / 1000), "k"), paste0(round(N / 1000), "k"))
    ) +
    scale_y_continuous(limits = c(0, 2.6), expand = c(0, 0)) +
    labs(x = "Rank", y = NULL) +
    theme_classic(base_size = 7.5) +
    theme(
      axis.text.y   = element_blank(),
      axis.ticks.y  = element_blank(),
      axis.line.y   = element_blank(),
      axis.text.x   = element_text(size = 5.5),
      axis.title.x  = element_text(size = 6, hjust = 0.5),
      plot.margin   = margin(0, 3, 3, 3)
    )

  p_es / p_bc + plot_layout(heights = c(3, 2))
}

# ─── Build 5 panels ──────────────────────────────────────────────────────────
cat("Generating panels ...\n")
MTORC1_GENES   <- c("GSK3B", "PDK1", "ADIPOR2", "ACLY", "GOT1", "GBE1")

panels <- lapply(SELECTED, function(pw) {
  cat(sprintf("  %s\n", pw))
  blue <- if (pw == "HALLMARK_MTORC1_SIGNALING") MTORC1_GENES else ADIPO_GENES
  make_panel(pw, blue_genes = blue)
})
names(panels) <- SELECTED

# ─── Arrange 3 panels in a row ───────────────────────────────────────────────
ordered_pws <- c(
  "GOBP_EXTERNAL_ENCAPSULATING_STRUCTURE_ORGANIZATION",
  "HALLMARK_ADIPOGENESIS",
  "HALLMARK_MTORC1_SIGNALING"
)

combined <- wrap_plots(panels[ordered_pws], nrow = 1) +
  plot_annotation(
    title    = paste0(CT, " — GSEA interaction term enrichment (supplementary)"),
    subtitle = paste0(
      "Ranked by \u0394logFC (obese \u2212 non-obese VAT/SAT); ",
      "NES / adj.p from limma interaction GSEA\n",
      "Orange ticks: fibrosis genes  |  Blue ticks: adipogenesis genes"
    ),
    theme = theme(
      plot.title    = element_text(size = 8, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 6.5, hjust = 0.5, color = "grey40")
    )
  )

# ─── Save ─────────────────────────────────────────────────────────────────────
out <- file.path(OUT_DIR,
  paste0(CT, "_GSEA_enrichment_5pathways_supp.pdf"))
ggsave(out, combined, width = 6, height = 3.2, device = "pdf")
cat(sprintf("\nSaved: %s (6 x 3.2 in)\n", basename(out)))
cat("[Done]\n")
