suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(edgeR)
  library(Matrix)
  library(Seurat)
  library(msigdbr)
})

# ─── Config ───────────────────────────────────────────────────────────────────
# Project root is 7 levels up from this script
args0 <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args0[grep("--file=", args0)])
project_root <- normalizePath(file.path(dirname(normalizePath(script_path)),
                                        paste(rep("..", 7), collapse = "/")))

OUT_DIR <- file.path(project_root, "data/35.vascular/04.SMC_module_score")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

METADATA_PATH   <- file.path(project_root, "data/08.integration_snRNA_object/01.freeze_version/metadata_after_libraryQC_after_ctQC_freeze.tsv")
SMC_COUNT_SAT   <- file.path(project_root, "data/35.vascular/SMC_counts/SAT")
SMC_COUNT_VAT   <- file.path(project_root, "data/35.vascular/SMC_counts/VAT")

GSEA_FILE <- file.path(project_root, "data/21.differential_expressed_genes/01.pseudobulk/04.only_VAT_vs_SAT_pseudobulk/linear_model/interation_term/tables/GSEA/SMC.GSEA.sig.GOBP_interaction.tsv")

ASIAN_RACES  <- c("Chinese", "Malay", "Indian")
MIN_CELLS    <- 10
COL_VAT      <- "#F3AF44"
COL_SAT      <- "#67B3AD"
COL_PLAIN    <- "#a6cee3"
DEPOT_COLORS <- c(SAT = COL_SAT, VAT = COL_VAT)
GROUP_COLORS <- c(
  "Non-obese"      = "#efedf5",
  "Obese"          = "#bcbddc",
  "Severely obese" = "#807dba",
  "Morbidly obese" = "#54278f"
)
BMI_LEVELS <- c("Non-obese", "Obese", "Severely obese", "Morbidly obese")

PATHWAY_S1 <- "smooth_muscle_contraction"
PATHWAY_S2 <- "collagen_fibril_organization"

# ─── Load full pathway gene sets from MSigDB ──────────────────────────────────
message("[1] Loading full GOBP gene sets from MSigDB (msigdbr)...")
msig_all <- msigdbr(species = "Homo sapiens")
get_msig <- function(gs) unique(msig_all$gene_symbol[msig_all$gs_name == gs])

SMOOTH_GENES   <- get_msig("GOBP_SMOOTH_MUSCLE_CONTRACTION")
COLLAGEN_GENES <- get_msig("GOBP_COLLAGEN_FIBRIL_ORGANIZATION")
message("  GOBP_SMOOTH_MUSCLE_CONTRACTION: ", length(SMOOTH_GENES), " genes")
message("  GOBP_COLLAGEN_FIBRIL_ORGANIZATION: ", length(COLLAGEN_GENES), " genes")

# Save gene lists
write.table(data.frame(gene = SMOOTH_GENES),
            file.path(OUT_DIR, "00a.GOBP_SMOOTH_MUSCLE_CONTRACTION_all_genes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(gene = COLLAGEN_GENES),
            file.path(OUT_DIR, "00b.GOBP_COLLAGEN_FIBRIL_ORGANIZATION_all_genes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

ALL_TARGET_GENES <- unique(c(SMOOTH_GENES, COLLAGEN_GENES))
message("  Total unique genes: ", length(ALL_TARGET_GENES))

GENESETS <- list(
  smooth_muscle_contraction    = SMOOTH_GENES,
  collagen_fibril_organization = COLLAGEN_GENES
)

# ─── Load metadata ────────────────────────────────────────────────────────────
message("[2] Loading metadata...")
meta <- fread(METADATA_PATH, sep = "\t", header = TRUE, data.table = FALSE)

cell_meta <- meta %>%
  select(barcode, donor_id, Depot, Race, BMI, Age, Gender, Diabetic) %>%
  filter(Race %in% ASIAN_RACES, !is.na(BMI))

donor_meta <- cell_meta %>%
  distinct(donor_id, .keep_all = TRUE) %>%
  mutate(
    T2DM          = ifelse(Diabetic == "T2DM", 1, 0),
    BMI_group     = case_when(
      BMI <  27.5 ~ "Non-obese",
      BMI <  32.5 ~ "Obese",
      BMI <  37.5 ~ "Severely obese",
      TRUE        ~ "Morbidly obese"
    ),
    BMI_group     = factor(BMI_group, levels = BMI_LEVELS),
    BMI_group_num = as.numeric(BMI_group)
  )

# ─── Load SMC counts (10x format) → Seurat → AddModuleScore ───────────────────
message("[3] Loading SMC counts and creating Seurat objects...")

load_smc_seurat <- function(count_dir, depot_name) {
  message("  ", depot_name, ": Reading from ", count_dir)
  mat <- Read10X(count_dir)
  seu <- CreateSeuratObject(counts = mat, min.cells = 0, min.features = 0)
  NormalizeData(seu, verbose = FALSE)
}

seu_sat <- load_smc_seurat(SMC_COUNT_SAT, "SAT")
seu_vat <- load_smc_seurat(SMC_COUNT_VAT, "VAT")

# Filter to Asian cells with valid BMI
valid_bc_sat <- cell_meta %>% filter(Depot == "SAT") %>% pull(barcode)
valid_bc_vat <- cell_meta %>% filter(Depot == "VAT") %>% pull(barcode)
seu_sat <- seu_sat[, intersect(colnames(seu_sat), valid_bc_sat)]
seu_vat <- seu_vat[, intersect(colnames(seu_vat), valid_bc_vat)]
message("  SAT cells: ", ncol(seu_sat), " | VAT cells: ", ncol(seu_vat))

# ─── Donor filtering ──────────────────────────────────────────────────────────
cnt_sat <- cell_meta %>%
  filter(Depot == "SAT", barcode %in% colnames(seu_sat)) %>%
  count(donor_id, name = "n_sat")
cnt_vat <- cell_meta %>%
  filter(Depot == "VAT", barcode %in% colnames(seu_vat)) %>%
  count(donor_id, name = "n_vat")
donors_keep <- inner_join(cnt_sat, cnt_vat, by = "donor_id") %>%
  filter(n_sat >= MIN_CELLS, n_vat >= MIN_CELLS) %>%
  pull(donor_id)
message("  Donors retained (SAT∩VAT >=", MIN_CELLS, " cells): ", length(donors_keep))

bc_sat_keep <- cell_meta %>% filter(Depot == "SAT", donor_id %in% donors_keep) %>% pull(barcode)
bc_vat_keep <- cell_meta %>% filter(Depot == "VAT", donor_id %in% donors_keep) %>% pull(barcode)
seu_sat <- seu_sat[, intersect(colnames(seu_sat), bc_sat_keep)]
seu_vat <- seu_vat[, intersect(colnames(seu_vat), bc_vat_keep)]

# ─── Seurat AddModuleScore ────────────────────────────────────────────────────
message("[4] Running AddModuleScore (Seurat)...")
seu_sat <- AddModuleScore(seu_sat, features = list(SMOOTH_GENES, COLLAGEN_GENES),
                          name = "smc_module", verbose = FALSE)
seu_vat <- AddModuleScore(seu_vat, features = list(SMOOTH_GENES, COLLAGEN_GENES),
                          name = "smc_module", verbose = FALSE)
# smc_module1 = smooth_muscle_contraction, smc_module2 = collagen_fibril_organization

# ─── Build per-cell table ─────────────────────────────────────────────────────
message("[5] Building per-cell score table...")

extract_scores <- function(seu, depot_name) {
  data.frame(
    cell_barcode                 = colnames(seu),
    Depot                        = depot_name,
    smooth_muscle_contraction    = seu$smc_module1,
    collagen_fibril_organization = seu$smc_module2
  ) %>%
    left_join(cell_meta %>% select(barcode, donor_id),
              by = c("cell_barcode" = "barcode")) %>%
    left_join(donor_meta %>% select(donor_id, BMI, BMI_group, BMI_group_num),
              by = "donor_id")
}

cell_scores <- bind_rows(
  extract_scores(seu_sat, "SAT"),
  extract_scores(seu_vat, "VAT")
) %>%
  mutate(Depot     = factor(Depot, levels = c("SAT", "VAT")),
         BMI_group = factor(BMI_group, levels = BMI_LEVELS))

# Save per-cell scores
fwrite(cell_scores,
       file.path(OUT_DIR, "02.SMC_module_score_per_cell.tsv.gz"),
       sep = "\t", quote = FALSE, row.names = FALSE, compress = "gzip")
message("  Saved: 02.SMC_module_score_per_cell.tsv.gz  (",
        nrow(cell_scores), " cells)")

# ─── Pseudobulk: mean module score per donor × depot ─────────────────────────
message("[6] Aggregating per donor...")
score_df <- cell_scores %>%
  group_by(donor_id, Depot) %>%
  summarise(
    smooth_muscle_contraction    = mean(smooth_muscle_contraction,    na.rm = TRUE),
    collagen_fibril_organization = mean(collagen_fibril_organization, na.rm = TRUE),
    n_cells = n(),
    .groups = "drop"
  ) %>%
  left_join(donor_meta %>%
              select(donor_id, BMI, BMI_group, BMI_group_num, Age, Gender, Race, T2DM),
            by = "donor_id") %>%
  filter(!is.na(BMI_group)) %>%
  mutate(Depot     = factor(Depot, levels = c("SAT", "VAT")),
         BMI_group = factor(BMI_group, levels = BMI_LEVELS))

write.table(score_df,
            file.path(OUT_DIR, "03.SMC_module_score_per_donor.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
message("  Saved: 03.SMC_module_score_per_donor.tsv")

# ─── LM + Spearman stats ─────────────────────────────────────────────────────
message("[7] Statistical tests...")

fmt_sci     <- function(p) formatC(p, format = "e", digits = 2)
score_long  <- score_df %>%
  pivot_longer(all_of(names(GENESETS)), names_to = "pathway", values_to = "score")

safe_lm_score <- function(df) {
  covars <- c("Age","Gender","T2DM","Race")
  keep <- covars[sapply(covars, function(v) {
    x <- df[[v]]
    if (is.numeric(x)) sd(x, na.rm = TRUE) > 0 else length(unique(na.omit(x))) >= 2
  })]
  fmla <- as.formula(paste("score ~ BMI_group_num",
                           if (length(keep) > 0) paste("+", paste(keep, collapse = " + ")) else ""))
  tryCatch({
    m  <- lm(fmla, data = df)
    ct <- summary(m)$coefficients
    sp <- cor.test(df$score, df$BMI_group_num, method = "spearman", exact = FALSE)
    data.frame(
      n_donors = nrow(df),
      lm_beta  = ct["BMI_group_num","Estimate"],
      lm_se    = ct["BMI_group_num","Std. Error"],
      lm_p     = ct["BMI_group_num","Pr(>|t|)"],
      sp_rho   = sp$estimate,
      sp_p     = sp$p.value
    )
  }, error = function(e)
    data.frame(n_donors = nrow(df), lm_beta = NA, lm_se = NA, lm_p = NA, sp_rho = NA, sp_p = NA)
  )
}

stats_res <- score_long %>%
  group_by(pathway, Depot) %>%
  group_modify(~ safe_lm_score(.x)) %>%
  ungroup() %>%
  mutate(
    lm_p_adj  = p.adjust(lm_p, method = "BH"),
    sp_p_adj  = p.adjust(sp_p, method = "BH"),
    sig_stars = case_when(
      is.na(lm_p)   ~ "?",
      lm_p < 0.001  ~ "***",
      lm_p < 0.01   ~ "**",
      lm_p < 0.05   ~ "*",
      lm_p < 0.1    ~ ".",
      TRUE          ~ "ns"
    ),
    label = paste0("beta=", round(lm_beta, 3), ", p=", fmt_sci(lm_p),
                   " (adj.p=", fmt_sci(lm_p_adj), ")",
                   "\nSpearman rho=", round(sp_rho, 3), ", p=", fmt_sci(sp_p))
  )

write.table(stats_res,
            file.path(OUT_DIR, "04.SMC_module_score_BMI_stats.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
message("\n  Stats:")
print(stats_res %>% select(pathway, Depot, n_donors, lm_beta, lm_p, sig_stars, sp_rho, sp_p) %>%
        as.data.frame())

n_labels <- score_long %>%
  distinct(pathway, donor_id, Depot, BMI_group) %>%
  count(pathway, Depot, BMI_group)

# ─── Boxplot helper ───────────────────────────────────────────────────────────
PATHWAY_LABELS <- c(
  smooth_muscle_contraction    = "GOBP SMOOTH MUSCLE\nCONTRACTION",
  collagen_fibril_organization = "GOBP COLLAGEN FIBRIL\nORGANIZATION"
)

make_boxplot <- function(pathway_name, fill_mode = c("bmi","plain")) {
  fill_mode <- match.arg(fill_mode)
  df   <- score_long  %>% filter(pathway == pathway_name)
  sdf  <- stats_res   %>% filter(pathway == pathway_name) %>%
    mutate(trend = ifelse(!is.na(lm_p) & lm_p < 0.05, "sig", "ns"))
  ndf  <- n_labels    %>% filter(pathway == pathway_name)
  df2  <- df %>% left_join(sdf %>% select(Depot, trend, lm_p), by = "Depot")

  p <- ggplot(df2, aes(x = BMI_group, y = score)) +
    {if (fill_mode == "bmi")
      geom_boxplot(aes(fill = BMI_group), outlier.size = 0.4,
                   width = 0.35, color = "grey30", linewidth = 0.3)
     else
      geom_boxplot(fill = COL_PLAIN, outlier.size = 0.4,
                   width = 0.35, color = "grey30", linewidth = 0.3)} +
    geom_jitter(width = 0.1, size = 0.3, color = "grey40", alpha = 0.5) +
    geom_smooth(data = df2 %>% filter(trend == "sig"),
                aes(x = BMI_group_num, y = score),
                method = "lm", formula = y ~ x,
                color = "red", fill = "red",
                linewidth = 0.5, alpha = 0.15, inherit.aes = FALSE) +
    geom_smooth(data = df2 %>% filter(trend == "ns"),
                aes(x = BMI_group_num, y = score),
                method = "lm", formula = y ~ x,
                color = "grey50", fill = "grey70", linetype = "dashed",
                linewidth = 0.5, alpha = 0.15, inherit.aes = FALSE) +
    geom_text(data = ndf,
              aes(x = BMI_group, y = -Inf, label = paste0("n=", n)),
              inherit.aes = FALSE, vjust = -0.3, size = 1.9, color = "grey40") +
    geom_text(data = sdf,
              aes(x = 2.5, y = Inf, label = label),
              inherit.aes = FALSE, vjust = 1.1, hjust = 0.5, size = 2.0, lineheight = 0.9) +
    facet_wrap(~ Depot, nrow = 1) +
    {if (fill_mode == "bmi")
      scale_fill_manual(values = GROUP_COLORS, name = "BMI group")
     else
      scale_fill_identity()} +
    labs(x = NULL, y = "Module score (Seurat AddModuleScore)",
         title = PATHWAY_LABELS[pathway_name]) +
    theme_classic(base_size = 9) +
    theme(
      axis.text.x      = element_text(size = 7, angle = 30, vjust = 1, hjust = 1),
      strip.background = element_blank(),
      strip.text       = element_text(face = "bold", size = 9),
      legend.position  = if (fill_mode == "bmi") "right" else "none",
      legend.key.size  = unit(0.3, "cm"),
      legend.text      = element_text(size = 6.5),
      legend.title     = element_text(size = 7, face = "bold"),
      plot.title       = element_text(size = 8.5, face = "bold")
    )
  p
}

# ─── Boxplots ─────────────────────────────────────────────────────────────────
message("[8] Plotting boxplots...")

plots_bmi   <- lapply(names(GENESETS), make_boxplot, fill_mode = "bmi")
plots_plain <- lapply(names(GENESETS), make_boxplot, fill_mode = "plain")

ggsave(file.path(OUT_DIR, "05a.SMC_module_score_boxplot_BMIcolor.pdf"),
       wrap_plots(plots_bmi, ncol = 1), width = 3.5, height = 4.5)
ggsave(file.path(OUT_DIR, "05b.SMC_module_score_boxplot_plaincolor.pdf"),
       wrap_plots(plots_plain, ncol = 1), width = 4, height = 4.5)
message("  Saved boxplots")

# ─── Correlation scatter plots ─────────────────────────────────────────────────
message("[9] Correlation scatter plots...")

# x = collagen (synthetic/fibrotic), y = smooth muscle contraction (contractile)
S1 <- names(GENESETS)[2]   # collagen_fibril_organization → x-axis
S2 <- names(GENESETS)[1]   # smooth_muscle_contraction    → y-axis
S1_lab <- "GOBP COLLAGEN FIBRIL ORGANIZATION\n(module score)"
S2_lab <- "GOBP SMOOTH MUSCLE CONTRACTION\n(module score)"

cor_stats <- score_df %>%
  group_by(Depot) %>%
  summarise(
    n     = n(),
    r     = cor(.data[[S1]], .data[[S2]], use = "complete.obs"),
    p_val = tryCatch(
      cor.test(.data[[S1]], .data[[S2]])$p.value,
      error = function(e) NA
    ),
    .groups = "drop"
  ) %>%
  mutate(
    p_adj = p.adjust(p_val, method = "BH"),
    label = paste0("Pearson R=", round(r, 3), "\np=", fmt_sci(p_val),
                   "\npadj=", fmt_sci(p_adj), ", n=", n)
  )

write.table(cor_stats,
            file.path(OUT_DIR, "06.SMC_module_score_correlation_stats.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# correlation by BMI group
cor_bygroup <- score_df %>%
  group_by(Depot, BMI_group) %>%
  summarise(
    n     = n(),
    r     = tryCatch(cor(.data[[S1]], .data[[S2]], use = "complete.obs"), error = function(e) NA),
    p_val = tryCatch(cor.test(.data[[S1]], .data[[S2]])$p.value, error = function(e) NA),
    .groups = "drop"
  )
write.table(cor_bygroup,
            file.path(OUT_DIR, "06b.SMC_correlation_by_depot_BMIgroup.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Label positions — upper right corner
label_pos <- score_df %>%
  group_by(Depot) %>%
  summarise(
    x = max(.data[[S1]], na.rm = TRUE),
    y = max(.data[[S2]], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(cor_stats %>% select(Depot, label), by = "Depot")

make_scatter <- function(color_mode = c("bmi","plain")) {
  color_mode <- match.arg(color_mode)
  p <- ggplot(score_df, aes(x = .data[[S1]], y = .data[[S2]])) +
    geom_smooth(method = "lm", formula = y ~ x,
                color = "grey45", fill = "grey80", linewidth = 0.6, alpha = 0.2) +
    geom_text(data = label_pos, aes(x = x, y = y, label = label),
              inherit.aes = FALSE, hjust = 1, vjust = 1.1,
              size = 2.2, color = "grey30", lineheight = 0.9) +
    facet_wrap(~ Depot, nrow = 1) +
    labs(x = S1_lab, y = S2_lab,
         title = "SMC module score correlation (by depot)") +
    theme_classic(base_size = 9) +
    theme(
      strip.background = element_blank(),
      strip.text       = element_text(face = "bold", size = 9),
      plot.title       = element_text(size = 9, face = "bold"),
      legend.key.size  = unit(0.3, "cm"),
      legend.text      = element_text(size = 6.5),
      legend.title     = element_text(size = 7, face = "bold")
    )

  if (color_mode == "bmi") {
    p + geom_point(aes(color = BMI_group), size = 1.5, alpha = 0.8) +
      scale_color_manual(values = GROUP_COLORS, name = "BMI group")
  } else {
    p + geom_point(color = COL_PLAIN, size = 1.5, alpha = 0.8)
  }
}

ggsave(file.path(OUT_DIR, "07a.SMC_module_score_correlation_BMIcolor.pdf"),
       make_scatter("bmi"), width = 5.5, height = 3)
ggsave(file.path(OUT_DIR, "07b.SMC_module_score_correlation_plaincolor.pdf"),
       make_scatter("plain"), width = 5, height = 3)
message("  Saved scatter plots")

# ─── Combined PDF ─────────────────────────────────────────────────────────────
pdf(file.path(OUT_DIR, "08.SMC_module_score_all_figures.pdf"), width = 6, height = 4)
for (p in c(plots_bmi, plots_plain)) print(p)
print(make_scatter("bmi"))
print(make_scatter("plain"))
dev.off()

message("\n=== Stats Summary ===")
print(stats_res %>%
        select(pathway, Depot, n_donors, lm_beta, lm_p, lm_p_adj, sig_stars) %>%
        as.data.frame())
message("\n=== Correlation ===")
print(cor_stats %>% as.data.frame())
message("\n[Done] Script 02 complete.")
