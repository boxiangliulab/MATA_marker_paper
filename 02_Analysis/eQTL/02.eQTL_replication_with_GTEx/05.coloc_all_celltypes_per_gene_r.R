#!/usr/bin/env Rscript
# 05.coloc_all_celltypes_per_gene_r.R
#
# For each cell type with coloc results, compute per-gene Pearson r and
# concordance between SAT and VAT using ALL intersection cis-SNPs.
# Categorise genes as Shared / SAT_specific / VAT_specific via coloc.
# Output: barplots of median per-gene r and concordance across cell types.

library(data.table)
library(ggplot2)

## =========================================================
## Paths
## =========================================================
COLOC_ROOT <- "/home/users/nus/e1101920/project/adipose/05_model/depot_eQTL/results/eQTL_SAT_VAT_coloc"
CIS_ROOT   <- "/data/projects/11003054/e1101943/project/01.adipose_sqtl/data/09.eQTL/03.Saige_QTL/02.output/03.formal_calling/ancestry_specific_version"
OUTPUT_DIR <- "/data/projects/11003054/e1101943/project/01.adipose_sqtl/data/09.eQTL/03.Saige_QTL/02.output/04.replication/coloc_based_correlation"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

## =========================================================
## Colors
## =========================================================
COL_SHARED <- "#1f78b4"
COL_SAT    <- "#67B3AD"
COL_VAT    <- "#F3AF44"

CT_COLORS <- c(
  "Adipocyte"              = "#56BF95",
  "Early_preadipocyte"     = "#c7e9c0",
  "Committed_preadipocyte" = "#a1d99b",
  "Areg"                   = "#CCCC4D",
  "IGFBP2+_cell"           = "#74c476",
  "Mesothelial"            = "#DAE220",
  "AEC"                    = "#F0B6DA",
  "CEC"                    = "#DE77AE",
  "VEC"                    = "#C51B7D",
  "LEC"                    = "#8E0052",
  "Pericyte"               = "#fcc5c0",
  "SMC2"                   = "#fa9fb5",
  "SMC"                    = "#f768a1",
  "SMC1"                   = "#f768a1",
  "Endometrial"            = "#883186",
  "cDC1"                   = "#a6bddb",
  "cDC2"                   = "#74a9cf",
  "Activate_DC"            = "#3690c0",
  "pDC"                    = "#0570b0",
  "Classical_monocyte"     = "#243494",
  "Non_classical_monocyte" = "#061C58",
  "Intermediate_monocyte"  = "#1C1894",
  "TIM4+_ATM"              = "#bcbddc",
  "PVM"                    = "#9e9ac8",
  "LAM"                    = "#807dba",
  "TIM4+_CD11c+_ATM"       = "#6a51a3",
  "B"                      = "#fdae61",
  "Plasmablast"            = "#993404",
  "Plasma_cell"            = "#f46d43",
  "Treg"                   = "#A7D5F7",
  "CD4+_T"                 = "#9BB9FB",
  "CD8+_T"                 = "#5087F3",
  "CD56dim_CD16+_NK"       = "#4157E4",
  "CD56+_NK"               = "#253CD1",
  "Mast"                   = "#44B3E1"
)

## =========================================================
## Helper: per-gene r and concordance from cis files
## =========================================================
per_gene_stats <- function(cis_root, celltype, gene) {
  sat_f <- file.path(cis_root, "SAT", celltype, paste0(gene, "_cis"))
  vat_f <- file.path(cis_root, "VAT", celltype, paste0(gene, "_cis"))
  if (!file.exists(sat_f) || !file.exists(vat_f)) return(NULL)

  sat_dt <- fread(sat_f, select = c("MarkerID", "BETA"), showProgress = FALSE)
  vat_dt <- fread(vat_f, select = c("MarkerID", "BETA"), showProgress = FALSE)

  merged <- merge(sat_dt, vat_dt, by = "MarkerID", suffixes = c("_SAT", "_VAT"))
  n <- nrow(merged)
  if (n < 3L) return(NULL)

  r   <- unname(cor(merged$BETA_SAT, merged$BETA_VAT))
  pct <- 100 * sum(sign(merged$BETA_SAT) == sign(merged$BETA_VAT)) / n
  list(r = r, pct_concordant = pct, n_snps = n)
}

## =========================================================
## Cell types to process
## =========================================================
celltypes <- list.dirs(COLOC_ROOT, full.names = FALSE, recursive = FALSE)
celltypes <- celltypes[celltypes != ""]
message("Cell types found: ", paste(celltypes, collapse = ", "))

## =========================================================
## Main loop
## =========================================================
all_gene_r <- list()

for (ct in celltypes) {
  message("\n=== ", ct, " ===")

  shared_f  <- file.path(COLOC_ROOT, ct, "colocalized.tsv")
  sat_f     <- file.path(COLOC_ROOT, ct, "SAT_specific.tsv")
  vat_f     <- file.path(COLOC_ROOT, ct, "VAT_specific.tsv")

  if (!file.exists(shared_f)) { message("  [SKIP] No coloc files."); next }

  shared_genes   <- fread(shared_f,  header = TRUE)$gene
  sat_spec_genes <- fread(sat_f,     header = TRUE)$gene
  vat_spec_genes <- fread(vat_f,     header = TRUE)$gene

  gene_cat <- rbindlist(list(
    data.table(gene = shared_genes,   category = "Shared"),
    data.table(gene = sat_spec_genes, category = "SAT_specific"),
    data.table(gene = vat_spec_genes, category = "VAT_specific")
  ))

  message("  Shared=", length(shared_genes),
          " SAT_specific=", length(sat_spec_genes),
          " VAT_specific=", length(vat_spec_genes))

  res <- rbindlist(lapply(seq_len(nrow(gene_cat)), function(i) {
    g    <- gene_cat$gene[i]
    cat_i <- gene_cat$category[i]
    s    <- per_gene_stats(CIS_ROOT, ct, g)
    if (is.null(s)) return(NULL)
    data.table(celltype = ct, gene = g, category = cat_i,
               r = s$r, pct_concordant = s$pct_concordant, n_snps = s$n_snps)
  }), fill = TRUE)

  if (!is.null(res) && nrow(res) > 0) {
    all_gene_r[[ct]] <- res
    message("  Done: ", nrow(res), " genes with data")
  }
}

## =========================================================
## Combine and save per-gene table
## =========================================================
all_dt <- rbindlist(all_gene_r, fill = TRUE)
all_dt[, category := factor(category,
                              levels = c("Shared", "SAT_specific", "VAT_specific"))]

fwrite(all_dt,
       file = file.path(OUTPUT_DIR, "all_celltypes_per_gene_r_concordance.tsv"),
       sep  = "\t")

## =========================================================
## Summary per cell type × category
## =========================================================
summary_dt <- all_dt[!is.na(r), .(
  mean_r        = mean(r),
  sd_r          = sd(r),
  mean_concordant = mean(pct_concordant),
  sd_conc       = sd(pct_concordant),
  n_genes       = .N
), by = .(celltype, category)]
summary_dt[, r_lo   := mean_r - sd_r]
summary_dt[, r_hi   := mean_r + sd_r]
summary_dt[, conc_lo := mean_concordant - sd_conc]
summary_dt[, conc_hi := mean_concordant + sd_conc]

summary_dt[, category := factor(category,
                                  levels = c("Shared", "SAT_specific", "VAT_specific"))]

fwrite(summary_dt,
       file = file.path(OUTPUT_DIR, "all_celltypes_summary_per_gene_r.tsv"),
       sep  = "\t")

## =========================================================
## Cell type order: by median Shared r descending
## =========================================================
ct_order <- summary_dt[category == "Shared"][order(-mean_r), celltype]
# add any ct missing from Shared
ct_order <- c(ct_order, setdiff(unique(summary_dt$celltype), ct_order))
summary_dt[, celltype := factor(celltype, levels = ct_order)]
all_dt[,     celltype := factor(celltype, levels = ct_order)]

## =========================================================
## Circle annotation data
## =========================================================
circle_dt <- data.table(
  celltype = factor(ct_order, levels = ct_order),
  color    = CT_COLORS[ct_order]
)
circle_dt[is.na(color), color := "grey50"]

## =========================================================
## Helper: barplot (mean ± SD) with circles below x-axis
## =========================================================
make_barplot <- function(data, y_col, q_lo, q_hi, y_lab, y_lo, y_hi,
                         title_str, circle_dt, ref_line = NULL) {

  step <- y_hi - y_lo
  cy   <- y_lo - step * 0.10
  ly   <- y_lo - step * 0.22

  ct_levels  <- levels(circle_dt$celltype)
  circle_ann <- data.table(
    celltype = factor(ct_levels, levels = ct_levels),
    y        = cy,  ly = ly,
    col      = circle_dt$color[match(ct_levels, circle_dt$celltype)]
  )

  fill_vals <- c("Shared" = COL_SHARED, "SAT_specific" = COL_SAT,
                 "VAT_specific" = COL_VAT,
                 setNames(circle_ann$col, as.character(circle_ann$celltype)))

  ggplot(data, aes(x = celltype, y = .data[[y_col]], fill = category)) +
    { if (!is.null(ref_line))
        geom_hline(yintercept = ref_line, linetype = "dashed",
                   color = "grey50", linewidth = 0.5) } +
    geom_col(position = position_dodge(0.8), width = 0.7) +
    geom_errorbar(aes(ymin = .data[[q_lo]], ymax = .data[[q_hi]]),
                  position = position_dodge(0.8), width = 0.25, linewidth = 0.4) +
    geom_point(data = circle_ann, aes(x = celltype, y = y, fill = celltype),
               inherit.aes = FALSE, shape = 21, size = 6, stroke = 0.6, color = "black") +
    geom_text(data = circle_ann, aes(x = celltype, y = ly, label = celltype),
              inherit.aes = FALSE, angle = 45, hjust = 1, size = 2.8, color = "black") +
    scale_fill_manual(values = fill_vals,
                      breaks = c("Shared", "SAT_specific", "VAT_specific")) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(ylim = c(y_lo, y_hi), clip = "off") +
    labs(x = NULL, y = y_lab, fill = "eQTL category", title = title_str) +
    theme_bw(base_size = 11) +
    theme(plot.title         = element_text(hjust = 0.5, face = "bold"),
          axis.text.x        = element_blank(),
          axis.ticks.x       = element_blank(),
          panel.grid.major.x = element_blank(),
          legend.position    = "top",
          plot.margin        = margin(5, 10, 80, 10))
}

## =========================================================
## Helper: boxplot with median highlighted and circles below x-axis
## =========================================================
make_boxplot <- function(data, y_col, y_lab, y_lo, y_hi,
                         title_str, circle_dt,
                         ref_line = NULL) {

  dw   <- 0.8          # dodge width
  bw   <- 0.65         # box width
  step <- y_hi - y_lo
  cy   <- y_lo - step * 0.10
  ly   <- y_lo - step * 0.22

  ct_levels  <- levels(circle_dt$celltype)
  circle_ann <- data.table(
    celltype = factor(ct_levels, levels = ct_levels),
    y        = cy,
    ly       = ly,
    col      = circle_dt$color[match(ct_levels, circle_dt$celltype)]
  )

  fill_vals <- c("Shared" = COL_SHARED, "SAT_specific" = COL_SAT,
                 "VAT_specific" = COL_VAT,
                 setNames(circle_ann$col, as.character(circle_ann$celltype)))

  ggplot(data, aes(x = celltype, y = .data[[y_col]], fill = category)) +
    { if (!is.null(ref_line))
        geom_hline(yintercept = ref_line, linetype = "dashed",
                   color = "grey50", linewidth = 0.5) } +
    geom_boxplot(position = position_dodge(dw), width = bw,
                 color = "grey30", linewidth = 0.35,
                 outlier.size = 0.4, outlier.alpha = 0.4,
                 fatten = NULL) +        # remove default median line
    stat_summary(aes(group = category),
                 fun = median, geom = "crossbar",
                 position = position_dodge(dw), width = bw * 0.95,
                 color = "white", linewidth = 1.0, fatten = 0) +
    geom_point(data = circle_ann, aes(x = celltype, y = y, fill = celltype),
               inherit.aes = FALSE, shape = 21, size = 6, stroke = 0.6, color = "black") +
    geom_text(data = circle_ann, aes(x = celltype, y = ly, label = celltype),
              inherit.aes = FALSE, angle = 45, hjust = 1, size = 2.8, color = "black") +
    scale_fill_manual(values = fill_vals,
                      breaks = c("Shared", "SAT_specific", "VAT_specific")) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(ylim = c(y_lo, y_hi), clip = "off") +
    labs(x = NULL, y = y_lab, fill = "eQTL category", title = title_str) +
    theme_bw(base_size = 11) +
    theme(plot.title         = element_text(hjust = 0.5, face = "bold"),
          axis.text.x        = element_blank(),
          axis.ticks.x       = element_blank(),
          panel.grid.major.x = element_blank(),
          legend.position    = "top",
          plot.margin        = margin(5, 10, 80, 10))
}

## =========================================================
## Plot 1: Pearson r boxplot
## =========================================================
p_r <- make_boxplot(
  data      = all_dt[!is.na(r)],
  y_col     = "r",
  y_lab     = "Per-gene Pearson r (SAT vs VAT, all cis-SNPs)",
  y_lo      = -0.5,
  y_hi      = 1.0,
  title_str = "SAT-VAT eQTL effect concordance by cell type (coloc categories)",
  circle_dt = circle_dt
)

ggsave(
  filename = file.path(OUTPUT_DIR, "boxplot_pearson_r_all_celltypes.pdf"),
  plot     = p_r,
  width    = max(7, length(ct_order) * 0.45),
  height   = 4.5,
  useDingbats = FALSE
)

## =========================================================
## Plot 2: Concordance boxplot
## =========================================================
p_c <- make_boxplot(
  data      = all_dt[!is.na(pct_concordant)],
  y_col     = "pct_concordant",
  y_lab     = "Per-gene concordance (% cis-SNPs, same direction)",
  y_lo      = 0,
  y_hi      = 100,
  title_str = "SAT-VAT eQTL direction concordance by cell type (coloc categories)",
  circle_dt = circle_dt,
  ref_line  = 50
)

ggsave(
  filename = file.path(OUTPUT_DIR, "boxplot_concordance_all_celltypes.pdf"),
  plot     = p_c,
  width    = max(7, length(ct_order) * 0.45),
  height   = 4.5,
  useDingbats = FALSE
)

## =========================================================
## Plot 3: Pearson r barplot (mean ± SD)
## =========================================================
p_r_bar <- make_barplot(
  data      = copy(summary_dt),
  y_col     = "mean_r",
  q_lo      = "r_lo",
  q_hi      = "r_hi",
  y_lab     = "Mean per-gene Pearson r \u00b1 SD (SAT vs VAT, all cis-SNPs)",
  y_lo      = -0.3,
  y_hi      = 1.0,
  title_str = "SAT-VAT eQTL effect concordance by cell type (coloc categories)",
  circle_dt = circle_dt
)

ggsave(
  filename = file.path(OUTPUT_DIR, "barplot_pearson_r_all_celltypes.pdf"),
  plot     = p_r_bar,
  width    = max(7, length(ct_order) * 0.45),
  height   = 4.5,
  useDingbats = FALSE
)

## =========================================================
## Plot 4: Concordance barplot (mean ± SD)
## =========================================================
p_c_bar <- make_barplot(
  data      = copy(summary_dt),
  y_col     = "mean_concordant",
  q_lo      = "conc_lo",
  q_hi      = "conc_hi",
  y_lab     = "Mean per-gene concordance \u00b1 SD (% cis-SNPs, same direction)",
  y_lo      = 30,
  y_hi      = 100,
  title_str = "SAT-VAT eQTL direction concordance by cell type (coloc categories)",
  circle_dt = circle_dt,
  ref_line  = 50
)

ggsave(
  filename = file.path(OUTPUT_DIR, "barplot_concordance_all_celltypes.pdf"),
  plot     = p_c_bar,
  width    = max(7, length(ct_order) * 0.45),
  height   = 4.5,
  useDingbats = FALSE
)

message("\nDone. Output in: ", OUTPUT_DIR)
