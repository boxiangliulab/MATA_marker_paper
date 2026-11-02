# 05.SMC_umap_plot.R
# Redraw SMC UMAP feature plots in R using rasterised points
# Output: 600 DPI PDF, small file size

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ggrastr)
  library(patchwork)
})

# ─── Paths ────────────────────────────────────────────────────────────────────
# Project root is 7 levels up from this script
args0 <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args0[grep("--file=", args0)])
project_root <- normalizePath(file.path(dirname(normalizePath(script_path)),
                                        paste(rep("..", 7), collapse = "/")))

UMAP_TSV  <- file.path(project_root, "data/35.vascular/05.SMC_umap/01.SMC_umap_coordinates_per_depot.tsv")
SCORE_TSV <- file.path(project_root, "data/35.vascular/04.SMC_module_score/02.SMC_module_score_per_cell.tsv.gz")
OUT_DIR   <- file.path(project_root, "data/35.vascular/05.SMC_umap")

# ─── Config ───────────────────────────────────────────────────────────────────
MODULE_LABELS <- c(
  smooth_muscle_contraction    = "Smooth muscle contraction",
  collagen_fibril_organization = "Collagen fibril organisation"
)
DEPOT_LABELS <- c(SAT = "SAT (subcutaneous)", VAT = "VAT (visceral)")
PT_SIZE  <- 0.3
PT_ALPHA <- 0.8

# ─── Load ─────────────────────────────────────────────────────────────────────
# UMAP coordinates from Python script (first col = cell_barcode, unnamed)
umap_df <- as.data.frame(fread(UMAP_TSV))
colnames(umap_df)[1] <- "cell_barcode"
umap_df <- umap_df[, c("cell_barcode", "UMAP_1", "UMAP_2", "Depot")]

# Seurat AddModuleScore (full MSigDB pathway genes) from script 02
score_df <- as.data.frame(fread(SCORE_TSV,
  select = c("cell_barcode", "smooth_muscle_contraction", "collagen_fibril_organization")))

df <- merge(umap_df, score_df, by = "cell_barcode", all.x = TRUE)
df$Depot <- factor(df$Depot, levels = c("SAT", "VAT"))
cat(sprintf("Loaded %d cells (UMAP + Seurat scores)\n", nrow(df)))

# ─── Helper: one feature plot panel ───────────────────────────────────────────
feature_panel <- function(sub, score_col, depot, score_label) {
  # plot low scores first so high-expressing cells visible on top
  sub <- sub[order(sub[[score_col]]), ]

  # shared colour limits: 1st–99th percentile
  ql <- quantile(sub[[score_col]], c(0.01, 0.99), na.rm = TRUE)

  ggplot(sub, aes(x = UMAP_1, y = UMAP_2, color = .data[[score_col]])) +
    rasterise(
      geom_point(size = PT_SIZE, alpha = PT_ALPHA, stroke = 0, shape = 16),
      dpi = 600
    ) +
    scale_color_viridis_c(
      option = "viridis",
      limits = ql,
      oob    = scales::squish,
      name   = score_label,
      guide  = guide_colorbar(
        title.position = "top",
        barwidth  = unit(2.2, "cm"),
        barheight = unit(0.22, "cm"),
        ticks     = FALSE
      )
    ) +
    coord_fixed() +
    labs(x = "UMAP 1", y = "UMAP 2",
         title = DEPOT_LABELS[[depot]]) +
    theme_classic(base_size = 8) +
    theme(
      plot.title       = element_text(size = 8, face = "bold", hjust = 0.5),
      axis.text        = element_blank(),
      axis.ticks       = element_blank(),
      axis.line        = element_line(linewidth = 0.3),
      legend.position  = "bottom",
      legend.title     = element_text(size = 7),
      legend.text      = element_text(size = 6.5),
      legend.margin    = margin(0, 0, 0, 0),
      plot.margin      = margin(2, 2, 2, 2)
    )
}

# ─── Build panels: rows = module, cols = depot ────────────────────────────────
plots <- list()
for (mod in names(MODULE_LABELS)) {
  for (dep in c("SAT", "VAT")) {
    sub <- df[df$Depot == dep, ]
    plots[[paste(mod, dep, sep = "_")]] <- feature_panel(sub, mod, dep, MODULE_LABELS[[mod]])
  }
}

# ─── Combined figure: 2 rows × 2 cols ────────────────────────────────────────
fig <- (plots[[1]] | plots[[2]]) / (plots[[3]] | plots[[4]]) +
  plot_annotation(
    title = "SMC module scores on per-depot UMAP",
    theme = theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5))
  )

out_combined <- file.path(OUT_DIR, "05.SMC_umap_module_score_R.pdf")
ggsave(out_combined, fig,
       width = 4.8, height = 4.8,
       device = cairo_pdf,
       dpi = 600)
cat(sprintf("Saved: %s (%.0f KB)\n", basename(out_combined),
            file.size(out_combined) / 1024))

# ─── Individual PDFs per module ───────────────────────────────────────────────
for (mod in names(MODULE_LABELS)) {
  p <- (plots[[paste(mod, "SAT", sep="_")]] | plots[[paste(mod, "VAT", sep="_")]]) +
    plot_annotation(
      title = MODULE_LABELS[[mod]],
      theme = theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5))
    )
  safe <- gsub("_", "-", mod)
  out  <- file.path(OUT_DIR, sprintf("06.SMC_umap_%s_R.pdf", safe))
  ggsave(out, p, width = 4.8, height = 2.6, device = cairo_pdf, dpi = 600)
  cat(sprintf("Saved: %s (%.0f KB)\n", basename(out), file.size(out) / 1024))
}

# ─── Load metadata: nCount_RNA + BMI group per cell ──────────────────────────
META_TSV  <- file.path(project_root, "data/08.integration_snRNA_object/01.freeze_version/metadata_after_libraryQC_after_ctQC_freeze.tsv")
COUNT_DIR <- file.path(project_root, "data/08.integration_snRNA_object/scvi/03.count_matrix")

meta <- as.data.frame(fread(META_TSV,
                             select = c("barcode", "donor_id", "BMI", "nCount_RNA")))
meta$BMI_group <- cut(meta$BMI,
                      breaks = c(-Inf, 27.5, 32.5, 37.5, Inf),
                      labels = c("Non-obese", "Obese", "Severely obese", "Morbidly obese"),
                      right  = FALSE)
df <- merge(df, meta, by.x = "cell_barcode", by.y = "barcode", all.x = TRUE)
df$BMI_group <- factor(df$BMI_group,
                       levels = c("Non-obese", "Obese", "Severely obese", "Morbidly obese"))
cat(sprintf("After metadata merge: %d / %d cells have BMI_group\n",
            sum(!is.na(df$BMI_group)), nrow(df)))

# ─── Load per-cell gene expression: MYH11, ACTA2, COL1A1 ─────────────────────
GENES3 <- c("MYH11", "ACTA2", "COL1A1")
GENE_LABELS <- c(MYH11  = "MYH11 (contractile)",
                 ACTA2  = "ACTA2 (contractile)",
                 COL1A1 = "COL1A1 (collagen)")

gene_list <- lapply(c("SAT", "VAT"), function(dep) {
  path <- file.path(COUNT_DIR, dep, "01.raw_count", "SMC.raw_count.tsv")
  hdr  <- fread(path, nrows = 0)
  cols <- c(1L, which(colnames(hdr) %in% GENES3))
  cnt  <- as.data.frame(fread(path, select = cols))
  colnames(cnt)[1] <- "cell_barcode"
  cnt$Depot <- dep
  cnt
})
gene_df <- do.call(rbind, gene_list)
df <- merge(df, gene_df, by = c("cell_barcode", "Depot"), all.x = TRUE)

# log-normalize: log1p(count / nCount_RNA × 10,000)
for (g in GENES3) {
  gcol <- paste0(g, "_lognorm")
  df[[gcol]] <- log1p(df[[g]] / df$nCount_RNA * 1e4)
}
cat(sprintf("Gene expression columns added for: %s\n", paste(GENES3, collapse = ", ")))

# ─── Gene expression feature plots ────────────────────────────────────────────
gene_plots <- list()
for (g in GENES3) {
  gcol <- paste0(g, "_lognorm")
  for (dep in c("SAT", "VAT")) {
    sub <- df[df$Depot == dep & !is.na(df[[gcol]]), ]
    sub <- sub[order(sub[[gcol]]), ]
    gene_plots[[paste(g, dep, sep = "_")]] <- feature_panel(sub, gcol, dep, GENE_LABELS[[g]])
  }
}

# combined 3-gene figure (3 rows × 2 cols)
gene_fig <- (gene_plots[["MYH11_SAT"]]  | gene_plots[["MYH11_VAT"]]) /
            (gene_plots[["ACTA2_SAT"]]  | gene_plots[["ACTA2_VAT"]]) /
            (gene_plots[["COL1A1_SAT"]] | gene_plots[["COL1A1_VAT"]]) +
  plot_annotation(
    title = "SMC contractile marker expression on per-depot UMAP",
    theme = theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5))
  )
out_genes <- file.path(OUT_DIR, "07.SMC_umap_genes_R.pdf")
ggsave(out_genes, gene_fig, width = 4.8, height = 7.0, device = cairo_pdf, dpi = 600)
cat(sprintf("Saved: %s (%.0f KB)\n", basename(out_genes), file.size(out_genes) / 1024))

# individual PDFs per gene (1 row × 2 cols)
for (g in GENES3) {
  p <- (gene_plots[[paste(g, "SAT", sep = "_")]] |
        gene_plots[[paste(g, "VAT", sep = "_")]]) +
    plot_annotation(title = GENE_LABELS[[g]],
                    theme = theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5)))
  out <- file.path(OUT_DIR, sprintf("07.SMC_umap_%s_R.pdf", g))
  ggsave(out, p, width = 4.8, height = 2.6, device = cairo_pdf, dpi = 600)
  cat(sprintf("Saved: %s (%.0f KB)\n", basename(out), file.size(out) / 1024))
}

# ─── BMI group facet UMAP (smooth muscle contraction score) ──────────────────
bmi_panel <- function(sub_dep, dep) {
  ql <- quantile(sub_dep$smooth_muscle_contraction, c(0.01, 0.99), na.rm = TRUE)
  ggplot(sub_dep, aes(x = UMAP_1, y = UMAP_2)) +
    rasterise(
      geom_point(aes(color = smooth_muscle_contraction),
                 size = PT_SIZE * 0.8, alpha = PT_ALPHA, stroke = 0, shape = 16),
      dpi = 600
    ) +
    scale_color_viridis_c(
      option = "viridis", limits = ql, oob = scales::squish,
      name  = "Contraction score",
      guide = guide_colorbar(title.position = "top",
                             barwidth  = unit(1.8, "cm"),
                             barheight = unit(0.2, "cm"),
                             ticks     = FALSE)
    ) +
    facet_wrap(~ BMI_group, nrow = 1) +
    coord_fixed() +
    labs(x = "UMAP 1", y = "UMAP 2", title = DEPOT_LABELS[[dep]]) +
    theme_classic(base_size = 7) +
    theme(
      plot.title       = element_text(size = 8, face = "bold", hjust = 0.5),
      strip.text       = element_text(size = 6.5, face = "bold"),
      strip.background = element_blank(),
      axis.text        = element_blank(),
      axis.ticks       = element_blank(),
      axis.line        = element_line(linewidth = 0.3),
      legend.position  = "bottom",
      legend.title     = element_text(size = 6.5),
      legend.text      = element_text(size = 6),
      plot.margin      = margin(2, 2, 6, 2)
    )
}

bmi_sat <- bmi_panel(df[df$Depot == "SAT" & !is.na(df$BMI_group), ], "SAT")
bmi_vat <- bmi_panel(df[df$Depot == "VAT" & !is.na(df$BMI_group), ], "VAT")

bmi_fig <- bmi_sat / bmi_vat +
  plot_annotation(
    title = "SMC smooth muscle contraction score by BMI group",
    theme = theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5))
  )
out_bmi <- file.path(OUT_DIR, "08.SMC_umap_BMI_facet_contraction_R.pdf")
ggsave(out_bmi, bmi_fig, width = 7.5, height = 5.2, device = cairo_pdf, dpi = 600)
cat(sprintf("Saved: %s (%.0f KB)\n", basename(out_bmi), file.size(out_bmi) / 1024))

cat("\n[Done]\n")
