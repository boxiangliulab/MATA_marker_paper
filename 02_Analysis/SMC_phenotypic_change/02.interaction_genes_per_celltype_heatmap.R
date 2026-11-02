library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(grid)

# ============================================================
# Paths
# ============================================================
# Project root is 7 levels up from this script
args0 <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args0[grep("--file=", args0)])
project_root <- normalizePath(file.path(dirname(normalizePath(script_path)),
                                        paste(rep("..", 7), collapse = "/")))

interaction_dir <- file.path(project_root, "data/21.differential_expressed_genes/01.pseudobulk/04.only_VAT_vs_SAT_pseudobulk/linear_model/interation_term/tables")
non_obese_dir   <- file.path(project_root, "data/21.differential_expressed_genes/01.pseudobulk/04.only_VAT_vs_SAT_pseudobulk/BMIgrp_non_obese/01.DEG/01.table")
obese_dir       <- file.path(project_root, "data/21.differential_expressed_genes/01.pseudobulk/04.only_VAT_vs_SAT_pseudobulk/BMIgrp_obese/01.DEG/01.table")
out_dir         <- file.path(project_root, "data/21.differential_expressed_genes/01.pseudobulk/04.only_VAT_vs_SAT_pseudobulk/linear_model/interation_term/plots/summary")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

celltypes <- c("SMC", "CEC", "AEC", "VEC", "LEC", "Pericyte")

# Cell type colors
ct_colors <- c(
  "Adipocyte" = "#56BF95",
  "Early_preadipocyte" = "#c7e9c0",
  "Committed_preadipocyte" = "#a1d99b",
  "Areg" = "#CCCC4D",
  "IGFBP2+_cell" = "#74c476",
  "Mesothelial" = "#DAE220",
  "AEC" = "#F0B6DA",
  "CEC" = "#DE77AE",
  "VEC" = "#C51B7D",
  "LEC" = "#8E0052",
  "Pericyte" = "#fcc5c0",
  "SMC2" = "#fa9fb5",
  "SMC" = "#f768a1",
  "SMC1" = "#f768a1",
  "Endometrial" = "#883186",
  "cDC1" = "#a6bddb",
  "cDC2" = "#74a9cf",
  "Activate_DC" = "#3690c0",
  "pDC" = "#0570b0",
  "Classical_monocyte" = "#243494",
  "Non_classical_monocyte" = "#061C58",
  "Intermediate_monocyte" = "#1C1894",
  "TIM4+_ATM" = "#bcbddc",
  "PVM" = "#9e9ac8",
  "LAM" = "#807dba",
  "TIM4+_CD11c+_ATM" = "#6a51a3",
  "B" = "#fdae61",
  "Plasmablast" = "#993404",
  "Plasma_cell" = "#f46d43",
  "Treg" = "#A7D5F7",
  "CD4+_T" = "#9BB9FB",
  "CD8+_T" = "#5087F3",
  "CD56dim_CD16+_NK" = "#4157E4",
  "CD56+_NK" = "#253CD1",
  "Mast" = "#44B3E1"
)

# Heatmap color gradient (blue=low, red=high)
heatmap_colors <- c("#313695", "#4575b4", "#74add1", "#abd9e9", "#e0f3f8",
                    "#ffffbf", "#fee090", "#fdae61", "#f46d43", "#d73027", "#a50026")

# ============================================================
# 1. Read interaction term genes per cell type
# ============================================================
interaction_genes <- list()
for (ct in celltypes) {
  f <- file.path(interaction_dir, paste0(ct, ".DepotBMI.sig.limma.tsv"))
  if (file.exists(f)) {
    df <- read.delim(f, stringsAsFactors = FALSE)
    if (nrow(df) > 0) {
      interaction_genes[[ct]] <- as.character(df$gene)
    }
  }
}

# ============================================================
# 2. Read DEG logFC for non-obese and obese (all genes table)
# ============================================================
read_deg <- function(dir_path, suffix, celltypes) {
  res <- list()
  for (ct in celltypes) {
    f <- file.path(dir_path, paste0(ct, ".", suffix))
    if (file.exists(f)) {
      df <- read.delim(f, stringsAsFactors = FALSE)
      res[[ct]] <- df %>% select(gene, logFC)
    }
  }
  res
}

non_obese_deg <- read_deg(non_obese_dir, "non_obese.paired.edgeR.all.tsv", celltypes)
obese_deg     <- read_deg(obese_dir,     "obese.paired.edgeR.all.tsv",     celltypes)

# ============================================================
# 3. Build combined data for circle heatmap
# ============================================================
combined <- list()
for (ct in celltypes) {
  genes <- interaction_genes[[ct]]
  if (is.null(genes) || length(genes) == 0) next

  no_df <- non_obese_deg[[ct]]
  ob_df <- obese_deg[[ct]]
  if (is.null(no_df) || is.null(ob_df)) next

  df <- data.frame(gene = genes, celltype = ct, stringsAsFactors = FALSE) %>%
    left_join(no_df %>% rename(logFC_nonobese = logFC), by = "gene") %>%
    left_join(ob_df %>% rename(logFC_obese = logFC), by = "gene")
  combined[[ct]] <- df
}
combined_df <- bind_rows(combined)
combined_df <- combined_df %>% filter(!is.na(logFC_nonobese) & !is.na(logFC_obese))

# Order by cell type
combined_df$celltype <- factor(combined_df$celltype, levels = celltypes)
combined_df <- combined_df %>% arrange(celltype)

# ============================================================
# 4. Plot 1 — Circular heatmap (sectors = cell types)
# ============================================================

# Only include celltypes that have data in combined_df
celltypes_heatmap <- celltypes[celltypes %in% unique(as.character(combined_df$celltype))]

# Prepare per-sector gene lists and logFC vectors
ct_vec <- as.character(combined_df$celltype)

# Global color range
all_vals <- c(combined_df$logFC_nonobese, combined_df$logFC_obese)
val_range <- quantile(abs(all_vals), 0.96, na.rm = TRUE)
col_fun <- colorRamp2(seq(-val_range, val_range, length.out = 11), heatmap_colors)

# Build per-sector data
sector_data <- split(combined_df, combined_df$celltype)
# xlim: each sector spans 0 to n_genes
sector_xlim <- do.call(rbind, lapply(celltypes_heatmap, function(ct) {
  n <- nrow(sector_data[[ct]])
  c(0, n)
}))
rownames(sector_xlim) <- celltypes_heatmap

pdf_path1 <- file.path(out_dir, "circle_heatmap_depot_logFC_by_BMI_v2.pdf")
message("Writing: ", pdf_path1)
message("celltypes_heatmap: ", paste(celltypes_heatmap, collapse = ", "))
message("combined_df nrow: ", nrow(combined_df))
if (file.exists(pdf_path1)) file.remove(pdf_path1)
pdf(pdf_path1, width = 16, height = 16)

circos.clear()
circos.par(
  gap.degree = 6,
  cell.padding = c(0, 0, 0, 0),
  start.degree = 90,
  track.margin = c(0.005, 0.005)
)

circos.initialize(factors = factor(celltypes_heatmap, levels = celltypes_heatmap),
                  xlim = sector_xlim)

# Track 1 (outermost): Gene labels
circos.track(
  factors = factor(celltypes_heatmap, levels = celltypes_heatmap),
  ylim = c(0, 1),
  track.height = 0.25,
  bg.border = NA,
  panel.fun = function(x, y) {
    ct <- CELL_META$sector.index
    df <- sector_data[[ct]]
    n <- nrow(df)
    for (i in seq_len(n)) {
      circos.text(i - 0.5, 0, df$gene[i],
                  facing = "clockwise", niceFacing = TRUE,
                  cex = 1.25, adj = c(0, 0.5))
    }
  }
)

# Track 2 (outer heatmap ring): Obese logFC
circos.track(
  factors = factor(celltypes_heatmap, levels = celltypes_heatmap),
  ylim = c(0, 1),
  track.height = 0.1,
  bg.border = NA,
  panel.fun = function(x, y) {
    ct <- CELL_META$sector.index
    df <- sector_data[[ct]]
    n <- nrow(df)
    for (i in seq_len(n)) {
      circos.rect(i - 1, 0, i, 1,
                  col = col_fun(df$logFC_obese[i]),
                  border = "white", lwd = 0.3)
    }
  }
)

# Track 3 (inner heatmap ring): Non-obese logFC
circos.track(
  factors = factor(celltypes_heatmap, levels = celltypes_heatmap),
  ylim = c(0, 1),
  track.height = 0.1,
  bg.border = NA,
  panel.fun = function(x, y) {
    ct <- CELL_META$sector.index
    df <- sector_data[[ct]]
    n <- nrow(df)
    for (i in seq_len(n)) {
      circos.rect(i - 1, 0, i, 1,
                  col = col_fun(df$logFC_nonobese[i]),
                  border = "white", lwd = 0.3)
    }
  }
)

# Track 4 (innermost): Cell type annotation bar — continuous per sector
circos.track(
  factors = factor(celltypes_heatmap, levels = celltypes_heatmap),
  ylim = c(0, 1),
  track.height = 0.05,
  bg.border = NA,
  panel.fun = function(x, y) {
    ct <- CELL_META$sector.index
    xlim <- CELL_META$xlim
    circos.rect(xlim[1], 0, xlim[2], 1, col = ct_colors[ct], border = NA)
  }
)

# Legends
lgd_heatmap <- Legend(
  col_fun = col_fun,
  title = "logFC",
  title_position = "topleft",
  direction = "vertical"
)

lgd_ring <- Legend(
  labels = c("Inner ring: Non-obese", "Outer ring: Obese"),
  title = "Heatmap ring",
  type = "grid",
  legend_gp = gpar(fill = c("grey70", "grey40"))
)

lgd_ct <- Legend(
  labels = celltypes_heatmap,
  title = "Cell type",
  legend_gp = gpar(fill = ct_colors[celltypes_heatmap]),
  type = "grid"
)

draw(lgd_heatmap, x = unit(0.88, "npc"), y = unit(0.88, "npc"), just = c("left", "top"))
draw(lgd_ring,    x = unit(0.88, "npc"), y = unit(0.58, "npc"), just = c("left", "top"))
draw(lgd_ct,      x = unit(0.88, "npc"), y = unit(0.45, "npc"), just = c("left", "top"))

circos.clear()
dev.off()

cat("Plot 1 saved.\n")

# ============================================================
# 5. Plot 2 — Bar plot of sig/tested proportion by direction
# ============================================================
interaction_all <- list()
for (ct in celltypes) {
  f <- file.path(interaction_dir, paste0(ct, ".DepotBMI.sig.limma.tsv"))
  if (file.exists(f)) {
    df <- read.delim(f, stringsAsFactors = FALSE)
    if (nrow(df) > 0) {
      df$gene <- as.character(df$gene)
      interaction_all[[ct]] <- df
    }
  }
}
interaction_all_df <- bind_rows(interaction_all)

# Tested genes: union of non_obese and obese edgeR-filtered genes per cell type
n_tested <- sapply(celltypes, function(ct) {
  genes_no <- if (!is.null(non_obese_deg[[ct]])) non_obese_deg[[ct]]$gene else character(0)
  genes_ob <- if (!is.null(obese_deg[[ct]])) obese_deg[[ct]]$gene else character(0)
  length(unique(c(genes_no, genes_ob)))
})
names(n_tested) <- celltypes

# Order by total proportion (sig / tested) descending
ct_total <- interaction_all_df %>%
  group_by(Celltype) %>%
  summarise(total = n(), .groups = "drop") %>%
  mutate(n_tested = n_tested[Celltype],
         total_prop = total / n_tested) %>%
  arrange(desc(total_prop))
ct_order <- ct_total$Celltype

bar_data <- interaction_all_df %>%
  mutate(direction = ifelse(beta > 0, "Positive", "Negative")) %>%
  group_by(Celltype, direction) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(
    n_tested = n_tested[Celltype],
    prop = count / n_tested,
    Celltype = factor(Celltype, levels = ct_order),
    plot_prop = ifelse(direction == "Negative", -prop, prop)
  )

p <- ggplot(bar_data, aes(x = Celltype, y = plot_prop, fill = direction)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = count),
            vjust = ifelse(bar_data$direction == "Positive", -0.5, 1.3),
            size = 4) +
  scale_fill_manual(values = c("Positive" = "#d73027", "Negative" = "#4575b4"),
                    name = "BMI*Depot beta") +
  geom_hline(yintercept = 0, linewidth = 0.5) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  labs(x = "Cell type", y = "Proportion of interaction genes (sig / tested)") +
  scale_x_discrete(labels = function(x) gsub("_", " ", x)) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )

pdf_path2 <- file.path(out_dir, "barplot_interaction_beta_direction_v2.pdf")
if (file.exists(pdf_path2)) file.remove(pdf_path2)
ggsave(pdf_path2,
       p, width = 8, height = 6)

cat("Plot 2 saved.\n")

# ============================================================
# Write figure captions
# ============================================================
captions <- c(
  "Figure: circle_heatmap_depot_logFC_by_BMI_v2.pdf",
  "Circular heatmap of BMI-dependent VAT-SAT expression differences across cell types.",
  "Each sector represents one cell type. Genes displayed are those with a significant",
  "BMI x Depot interaction term from the limma linear model (FDR < 0.05).",
  "The outer heatmap ring shows logFC (VAT vs SAT) in obese donors;",
  "the inner heatmap ring shows logFC (VAT vs SAT) in non-obese donors.",
  "Color scale (blue-to-red): negative = lower in VAT, positive = higher in VAT.",
  "The innermost bar indicates cell type identity.",
  "Only cell types with both non-obese and obese edgeR expression data are shown.",
  "",
  "Figure: barplot_interaction_beta_direction_v2.pdf",
  "Bar plot showing the proportion of BMI-dependent depot interaction genes per cell type.",
  "Y-axis: proportion = number of significant interaction genes / number of tested genes.",
  "Tested genes are defined as the union of genes passing edgeR filtering in the",
  "non-obese and obese paired edgeR analyses (edgeR.all.tsv files).",
  "Bars are split by the sign of the DepotVAT x BMI interaction beta coefficient:",
  "positive beta (red, upward bars) = VAT-SAT difference increases with BMI;",
  "negative beta (blue, downward bars) = VAT-SAT difference decreases with BMI.",
  "Numbers on bars indicate the raw count of significant interaction genes.",
  "Cell types are ordered by total proportion (positive + negative) in descending order."
)
writeLines(captions, file.path(out_dir, "figure_captions.txt"))

cat("Done!\n")
