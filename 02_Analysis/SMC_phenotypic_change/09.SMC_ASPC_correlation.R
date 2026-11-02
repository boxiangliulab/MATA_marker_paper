suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(cowplot)
})

# ─── Config ───────────────────────────────────────────────────────────────────
# Project root is 7 levels up from this script
args0 <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args0[grep("--file=", args0)])
project_root <- normalizePath(file.path(dirname(normalizePath(script_path)),
                                        paste(rep("..", 7), collapse = "/")))

OUT_DIR <- file.path(project_root, "data/35.vascular/06.SMC_ASPC_correlation")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

SMC_CELL_FILE   <- file.path(project_root, "data/35.vascular/04.SMC_module_score/02.SMC_module_score_per_cell.tsv.gz")
SMC_DONOR_FILE  <- file.path(project_root, "data/35.vascular/04.SMC_module_score/03.SMC_module_score_per_donor.tsv")
ASPC_CELL_FILE  <- file.path(project_root, "data/30.adipogenesis/02.committed_aspc_module_score_per_cell.tsv.gz")

MIN_CELLS <- 10

GROUP_COLORS <- c(
  "Non-obese"      = "#efedf5",
  "Obese"          = "#bcbddc",
  "Severely obese" = "#807dba",
  "Morbidly obese" = "#54278f"
)
BMI_LEVELS   <- c("Non-obese", "Obese", "Severely obese", "Morbidly obese")
DEPOT_COLORS <- c(SAT = "#67B3AD", VAT = "#F3AF44")

FIBROSIS_MODULES   <- c("Eisinger2023", "NABA_Collagens", "Collagen_Formation")
ADIPOGENESIS_MODULES <- c("HALLMARK_ADIPOGENESIS", "REACTOME_ADIPOGENESIS", "GOBP_WFC", "KEGG_PPAR")
ALL_ASPC_MODULES   <- c(FIBROSIS_MODULES, ADIPOGENESIS_MODULES)

# ─── 1. SMC cell proportions per donor (SAT) ──────────────────────────────────
message("[1] Computing SMC contractile/collagen proportions (SAT)...")
smc_cells <- fread(SMC_CELL_FILE, sep = "\t", header = TRUE, data.table = FALSE)

smc_prop <- smc_cells %>%
  filter(Depot == "SAT") %>%
  mutate(dominant = ifelse(smooth_muscle_contraction >= collagen_fibril_organization,
                           "contractile", "collagen")) %>%
  group_by(donor_id) %>%
  summarise(
    n_smc            = n(),
    prop_contractile = mean(dominant == "contractile"),
    prop_collagen    = mean(dominant == "collagen"),
    .groups = "drop"
  ) %>%
  filter(n_smc >= MIN_CELLS)

# Grab BMI_group from SMC cell file
bmi_lut <- smc_cells %>%
  distinct(donor_id, BMI, BMI_group) %>%
  mutate(BMI_group = factor(BMI_group, levels = BMI_LEVELS))

smc_prop <- smc_prop %>% left_join(bmi_lut, by = "donor_id")
message("  SMC donors (SAT, >=", MIN_CELLS, " cells): ", nrow(smc_prop))

# ─── 2. ASPC module scores per donor (SAT + VAT) ──────────────────────────────
message("[2] Aggregating ASPC module scores per donor (both depots)...")
aspc_cells <- fread(ASPC_CELL_FILE, sep = "\t", header = TRUE, data.table = FALSE)

aspc_donor_all <- aspc_cells %>%
  group_by(donor_id, Depot) %>%
  filter(n() >= MIN_CELLS) %>%
  summarise(
    across(all_of(ALL_ASPC_MODULES), ~ mean(.x, na.rm = TRUE)),
    n_aspc = n(),
    .groups = "drop"
  )

aspc_donor <- aspc_donor_all %>% filter(Depot == "SAT")
message("  ASPC donors (SAT, >=", MIN_CELLS, " cells): ", nrow(aspc_donor))
message("  ASPC donors (VAT, >=", MIN_CELLS, " cells): ",
        nrow(aspc_donor_all %>% filter(Depot == "VAT")))

# ─── 3. Merge ─────────────────────────────────────────────────────────────────
message("[3] Merging...")
merged <- inner_join(smc_prop, aspc_donor, by = "donor_id")
message("  Final donors: ", nrow(merged))

write.table(merged, file.path(OUT_DIR, "01.SMC_ASPC_correlation_data.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ─── 3b. SMC module scores per donor (both depots) ────────────────────────────
message("[3b] Loading SMC per-donor module scores (both depots)...")
smc_donor_all <- fread(SMC_DONOR_FILE, sep = "\t", header = TRUE, data.table = FALSE) %>%
  filter(n_cells >= MIN_CELLS) %>%
  select(donor_id, Depot, smooth_muscle_contraction, collagen_fibril_organization)
message("  SMC donors SAT: ", sum(smc_donor_all$Depot == "SAT"),
        " | VAT: ", sum(smc_donor_all$Depot == "VAT"))

smc_donor <- smc_donor_all %>% filter(Depot == "SAT") %>% select(-Depot)

merged_score <- inner_join(smc_donor, aspc_donor, by = "donor_id") %>%
  left_join(bmi_lut, by = "donor_id") %>% mutate(Depot = "SAT")

merged_score_vat <- inner_join(
  smc_donor_all %>% filter(Depot == "VAT") %>% select(-Depot),
  aspc_donor_all %>% filter(Depot == "VAT") %>% select(-Depot),
  by = "donor_id"
) %>% left_join(bmi_lut, by = "donor_id") %>% mutate(Depot = "VAT")

message("  Merged SAT: ", nrow(merged_score), " | VAT: ", nrow(merged_score_vat), " donors")

write.table(bind_rows(merged_score, merged_score_vat),
            file.path(OUT_DIR, "01b.SMC_score_ASPC_correlation_data.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ─── 4. Correlation statistics ────────────────────────────────────────────────
message("[4] Computing correlations...")

smc_vars  <- c("prop_contractile", "prop_collagen")
cor_stats <- expand.grid(smc_var = smc_vars, aspc_module = ALL_ASPC_MODULES,
                         stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(
    r     = cor.test(merged[[smc_var]], merged[[aspc_module]], method = "pearson")$estimate,
    p     = cor.test(merged[[smc_var]], merged[[aspc_module]], method = "pearson")$p.value,
    n     = sum(!is.na(merged[[smc_var]]) & !is.na(merged[[aspc_module]]))
  ) %>%
  ungroup() %>%
  mutate(
    padj       = p.adjust(p, method = "BH"),
    sig        = case_when(padj < 0.001 ~ "***", padj < 0.01 ~ "**",
                           padj < 0.05 ~ "*", TRUE ~ ""),
    module_type = ifelse(aspc_module %in% FIBROSIS_MODULES, "Fibrosis", "Adipogenesis")
  )

write.table(cor_stats, file.path(OUT_DIR, "02a.SMC_prop_ASPC_correlation_stats.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

message("\n  [Proportion] Correlation results:")
print(cor_stats %>% select(smc_var, aspc_module, r, p, padj, sig) %>% as.data.frame())

# ── Module score correlations ──
smc_score_vars <- c("smooth_muscle_contraction", "collagen_fibril_organization")
cor_stats_score <- expand.grid(smc_var = smc_score_vars, aspc_module = ALL_ASPC_MODULES,
                               stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(
    r = cor.test(merged_score[[smc_var]], merged_score[[aspc_module]], method = "pearson")$estimate,
    p = cor.test(merged_score[[smc_var]], merged_score[[aspc_module]], method = "pearson")$p.value,
    n = sum(!is.na(merged_score[[smc_var]]) & !is.na(merged_score[[aspc_module]]))
  ) %>%
  ungroup() %>%
  mutate(
    padj       = p.adjust(p, method = "BH"),
    sig        = case_when(padj < 0.001 ~ "***", padj < 0.01 ~ "**",
                           padj < 0.05 ~ "*", TRUE ~ ""),
    module_type = ifelse(aspc_module %in% FIBROSIS_MODULES, "Fibrosis", "Adipogenesis")
  )

write.table(cor_stats_score, file.path(OUT_DIR, "02b.SMC_score_ASPC_correlation_stats.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

message("\n  [Module score] Correlation results:")
print(cor_stats_score %>% select(smc_var, aspc_module, r, p, padj, sig) %>% as.data.frame())

# ── Also compute padj for VAT ──
cor_stats_score_vat <- expand.grid(smc_var = smc_score_vars, aspc_module = ALL_ASPC_MODULES,
                                    stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(
    r = cor.test(merged_score_vat[[smc_var]], merged_score_vat[[aspc_module]], method = "pearson")$estimate,
    p = cor.test(merged_score_vat[[smc_var]], merged_score_vat[[aspc_module]], method = "pearson")$p.value,
    n = sum(!is.na(merged_score_vat[[smc_var]]) & !is.na(merged_score_vat[[aspc_module]]))
  ) %>%
  ungroup() %>%
  mutate(
    padj = p.adjust(p, method = "BH"),
    sig  = case_when(padj < 0.001 ~ "***", padj < 0.01 ~ "**", padj < 0.05 ~ "*", TRUE ~ ""),
    module_type = ifelse(aspc_module %in% FIBROSIS_MODULES, "Fibrosis", "Adipogenesis")
  )

write.table(cor_stats_score_vat, file.path(OUT_DIR, "02c.SMC_score_ASPC_correlation_stats_VAT.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
message("\n  [Module score] VAT correlation results:")
print(cor_stats_score_vat %>% select(smc_var, aspc_module, r, p, padj, sig) %>% as.data.frame())

# Build padj lookup: list by depot → named vector [smc_var:aspc_module] = padj
make_padj_lut <- function(df) {
  setNames(df$padj, paste(df$smc_var, df$aspc_module, sep = ":"))
}
padj_lut <- list(SAT = make_padj_lut(cor_stats_score),
                 VAT = make_padj_lut(cor_stats_score_vat))

# ─── 5. Heatmap of correlations ───────────────────────────────────────────────
message("[5] Plotting correlation heatmap (module scores, SAT + VAT)...")

make_heat <- function(df, depot_label) {
  df <- df %>%
    mutate(
      aspc_module = factor(aspc_module, levels = rev(ALL_ASPC_MODULES)),
      smc_var_lab = recode(smc_var,
                           smooth_muscle_contraction    = "Contraction\nscore",
                           collagen_fibril_organization = "Collagen\nscore")
    )
  r_lim <- max(abs(df$r)) * 1.05
  ggplot(df, aes(x = smc_var_lab, y = aspc_module, fill = r)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = ifelse(sig != "", sprintf("%.2f%s", r, sig),
                                 sprintf("%.2f", r))),
              size = 2.2, color = "black") +
    scale_fill_gradientn(
      colors = c("#543005","#bf812d","#f6e8c3","#f5f5f5","#c7eae5","#35978f","#003c30"),
      limits = c(-r_lim, r_lim), name = "Pearson r"
    ) +
    facet_grid(module_type ~ ., scales = "free_y", space = "free_y") +
    labs(x = NULL, y = NULL, title = depot_label) +
    theme_classic(base_size = 7) +
    theme(
      plot.title       = element_text(size = 8, face = "bold", hjust = 0.5),
      strip.background = element_blank(),
      strip.text       = element_text(size = 6.5, face = "bold"),
      legend.position  = "right",
      legend.key.height = unit(20, "pt"),
      legend.key.width  = unit(5, "pt"),
      legend.text      = element_text(size = 5.5),
      legend.title     = element_text(size = 6),
      axis.text        = element_text(size = 6),
      axis.line        = element_blank(),
      axis.ticks       = element_blank()
    )
}

# Shared colour scale across SAT and VAT
r_lim_global <- max(abs(c(cor_stats_score$r, cor_stats_score_vat$r))) * 1.05

heat_sat <- make_heat(cor_stats_score,     "SAT") +
  scale_fill_gradientn(
    colors = c("#543005","#bf812d","#f6e8c3","#f5f5f5","#c7eae5","#35978f","#003c30"),
    limits = c(-r_lim_global, r_lim_global), name = "Pearson r"
  )
heat_vat <- make_heat(cor_stats_score_vat, "VAT") +
  scale_fill_gradientn(
    colors = c("#543005","#bf812d","#f6e8c3","#f5f5f5","#c7eae5","#35978f","#003c30"),
    limits = c(-r_lim_global, r_lim_global), name = "Pearson r"
  )

ggsave(file.path(OUT_DIR, "03.SMC_ASPC_correlation_heatmap.pdf"),
       heat_sat + heat_vat + plot_layout(guides = "collect"),
       width = 5.5, height = 3.2, device = cairo_pdf)
message("  Saved: 03.SMC_ASPC_correlation_heatmap.pdf (5.5 x 3.2 in)")

# ─── 6. Scatter plots: SMC module score (x) vs ASPC module score (y) ──────────
# Corresponding pairs: collagen SMC ~ fibrosis ASPC | contraction SMC ~ adipogenesis ASPC
# SAT and VAT plotted separately
message("[6] Plotting scatter plots (module score, corresponding pairs, SAT & VAT separate)...")

clean_lab <- function(x) gsub("_", " ", x)

make_scatter_score <- function(dat, x_col, y_col, x_lab, y_lab, pt_color,
                               base_sz = 11, ann_sz = 3.5, padj_val = NULL) {
  cr   <- cor.test(dat[[x_col]], dat[[y_col]], method = "pearson")
  lbl  <- sprintf("r = %.2f\np = %.2e", cr$estimate, cr$p.value)
  if (!is.null(padj_val))
    lbl <- sprintf("r = %.2f\np = %.2e\nFDR = %.2e", cr$estimate, cr$p.value, padj_val)
  ggplot(dat, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    geom_smooth(method = "lm", se = TRUE, color = "grey50",
                linewidth = 0.5, linetype = "dashed", alpha = 0.15) +
    geom_point(size = 1.5, alpha = 0.8, color = pt_color) +
    annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.4,
             label = lbl, size = ann_sz, color = "black", lineheight = 1.1) +
    labs(x = x_lab, y = y_lab) +
    theme_classic(base_size = base_sz) +
    theme(axis.title = element_text(size = base_sz - 1),
          axis.text  = element_text(size = base_sz - 2))
}

make_depot_fig <- function(dat, pt_color, depot_label, lut) {
  row1_plots <- lapply(FIBROSIS_MODULES, function(m)
    make_scatter_score(dat, "collagen_fibril_organization", m,
                       "SMC collagen score", clean_lab(m), pt_color,
                       base_sz = 11, ann_sz = 3.5,
                       padj_val = lut[[paste0("collagen_fibril_organization:", m)]]))
  row2_plots <- lapply(ADIPOGENESIS_MODULES, function(m)
    make_scatter_score(dat, "smooth_muscle_contraction", m,
                       "SMC contraction score", clean_lab(m), pt_color,
                       base_sz = 11, ann_sz = 3.5,
                       padj_val = lut[[paste0("smooth_muscle_contraction:", m)]]))

  # Row 1: 3 fibrosis panels (each 4x2.5) → total row width 12
  row1 <- wrap_plots(row1_plots, nrow = 1) +
    plot_annotation(title = paste0("SMC Collagen ~ ASPC Fibrosis (", depot_label, ")"),
                    theme = theme(plot.title = element_text(size = 9, face = "bold")))
  # Row 2: 4 adipogenesis panels (same total width 12, each ~3x2.5)
  row2 <- wrap_plots(row2_plots, nrow = 1) +
    plot_annotation(title = paste0("SMC Contraction ~ ASPC Adipogenesis (", depot_label, ")"),
                    theme = theme(plot.title = element_text(size = 9, face = "bold")))
  row1 / row2
}

fig_sat <- make_depot_fig(merged_score,     DEPOT_COLORS["SAT"], "SAT", padj_lut[["SAT"]])
fig_vat <- make_depot_fig(merged_score_vat, DEPOT_COLORS["VAT"], "VAT", padj_lut[["VAT"]])

# Width = 3 panels × 4 in = 12; Height = 2 rows × 2.5 in + 0.5 titles = 5.5
ggsave(file.path(OUT_DIR, "04a.SMC_ASPC_scatter_SAT.pdf"),
       fig_sat, width = 12, height = 3.7, device = cairo_pdf)
ggsave(file.path(OUT_DIR, "04b.SMC_ASPC_scatter_VAT.pdf"),
       fig_vat, width = 12, height = 3.7, device = cairo_pdf)
message("  Saved: 04a/04b scatter PDFs (12 x 3.7 in)")

# ─── 7. Collagen ~ Adipogenesis scatter (2×2 square panels) ───────────────────
message("[7] Plotting collagen ~ adipogenesis scatter (square panels)...")

make_collagen_adipo_fig <- function(dat, pt_color, depot_label, lut) {
  plots <- lapply(ADIPOGENESIS_MODULES, function(m)
    make_scatter_score(dat, "collagen_fibril_organization", m,
                       "SMC collagen score", clean_lab(m), pt_color,
                       base_sz = 10, ann_sz = 3.2,
                       padj_val = lut[[paste0("collagen_fibril_organization:", m)]]))
  wrap_plots(plots, nrow = 2, ncol = 2) +
    plot_annotation(
      title = paste0("SMC Collagen ~ ASPC Adipogenesis (", depot_label, ")"),
      theme = theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5))
    )
}

fig_col_adipo_sat <- make_collagen_adipo_fig(merged_score,     DEPOT_COLORS["SAT"], "SAT", padj_lut[["SAT"]])
fig_col_adipo_vat <- make_collagen_adipo_fig(merged_score_vat, DEPOT_COLORS["VAT"], "VAT", padj_lut[["VAT"]])

# 2×2 grid: 2 panels × 3 in wide, 2 rows × 3 in tall = 6×6 square
ggsave(file.path(OUT_DIR, "05a.SMC_collagen_adipo_scatter_SAT.pdf"),
       fig_col_adipo_sat, width = 6, height = 6, device = cairo_pdf)
ggsave(file.path(OUT_DIR, "05b.SMC_collagen_adipo_scatter_VAT.pdf"),
       fig_col_adipo_vat, width = 6, height = 6, device = cairo_pdf)
message("  Saved: 05a/05b collagen~adipogenesis PDFs (6 x 6 in)")
message("[Done]")
