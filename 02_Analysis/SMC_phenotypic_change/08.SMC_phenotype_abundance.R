# 05.SMC_phenotype_abundance.R
# Per-donor proportion of Contractile vs Collagen/ECM SMC cells by BMI group.
# Classification: per-cell phenotype = whichever module score is higher.
# LM: proportion ~ BMI_group_num (+ covariates) per phenotype per depot.
# Output: two PDFs (one per depot), each with two boxplot panels (one per phenotype).

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
})

# ─── Paths ────────────────────────────────────────────────────────────────────
# Per-cell Seurat AddModuleScore (full MSigDB pathway genes) from script 02
# Project root is 7 levels up from this script
args0 <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args0[grep("--file=", args0)])
project_root <- normalizePath(file.path(dirname(normalizePath(script_path)),
                                        paste(rep("..", 7), collapse = "/")))

SCORE_TSV  <- file.path(project_root, "data/35.vascular/04.SMC_module_score/02.SMC_module_score_per_cell.tsv.gz")
META_TSV   <- file.path(project_root, "data/08.integration_snRNA_object/01.freeze_version/metadata_after_libraryQC_after_ctQC_freeze.tsv")
OUT_DIR    <- file.path(project_root, "data/35.vascular/05.SMC_umap")

BMI_LEVELS  <- c("Non-obese", "Obese", "Severely obese", "Morbidly obese")
DEPOT_COLORS <- c(SAT = "#67B3AD", VAT = "#F3AF44")
PHENO_COLORS <- c(Contractile = "#2171b5", "Collagen/ECM" = "#e6550d")
COVARS       <- c("Age", "Gender", "T2DM", "Race")

# ─── Load data ────────────────────────────────────────────────────────────────
# Per-cell Seurat scores (consistent with script 02 BMI LM and ASPC correlation)
df <- as.data.frame(fread(SCORE_TSV,
  select = c("cell_barcode", "Depot", "donor_id",
             "smooth_muscle_contraction", "collagen_fibril_organization")))

# BMI + covariates from metadata
meta <- as.data.frame(fread(META_TSV,
  select = c("barcode", "BMI", "Age", "Gender", "Race", "Diabetic")))
meta <- meta %>%
  mutate(
    T2DM      = ifelse(Diabetic == "T2DM", 1L, 0L),
    BMI_group = cut(BMI, breaks = c(-Inf, 27.5, 32.5, 37.5, Inf),
                    labels = BMI_LEVELS, right = FALSE),
    BMI_group = factor(BMI_group, levels = BMI_LEVELS),
    BMI_group_num = as.integer(BMI_group)
  )

df <- merge(df, meta, by.x = "cell_barcode", by.y = "barcode", all.x = TRUE)

# ─── Classify cells ───────────────────────────────────────────────────────────
df$phenotype <- ifelse(
  df$smooth_muscle_contraction > df$collagen_fibril_organization,
  "Contractile", "Collagen/ECM"
)
df$phenotype <- factor(df$phenotype, levels = c("Contractile", "Collagen/ECM"))

# ─── Per-donor proportions ────────────────────────────────────────────────────
prop_raw <- df %>%
  filter(!is.na(donor_id), !is.na(BMI_group)) %>%
  group_by(donor_id, Depot, BMI_group, BMI_group_num,
           Age, Gender, T2DM, Race) %>%
  summarise(
    n_total       = n(),
    n_contractile = sum(phenotype == "Contractile"),
    n_collagen    = sum(phenotype == "Collagen/ECM"),
    .groups = "drop"
  ) %>%
  mutate(
    prop_contractile = n_contractile / n_total,
    prop_collagen    = n_collagen    / n_total
  )

# long format for plotting
prop_long <- prop_raw %>%
  tidyr::pivot_longer(
    cols      = c(prop_contractile, prop_collagen),
    names_to  = "phenotype",
    values_to = "proportion"
  ) %>%
  mutate(phenotype = recode(phenotype,
    prop_contractile = "Contractile",
    prop_collagen    = "Collagen/ECM"
  ),
  phenotype = factor(phenotype, levels = c("Contractile", "Collagen/ECM")))

cat(sprintf("Donors: %d, BMI groups: %s\n",
            length(unique(prop_long$donor_id)),
            paste(table(unique(prop_raw[, c("donor_id","BMI_group")])$BMI_group),
                  collapse = " / ")))

# ─── LM per phenotype per depot ───────────────────────────────────────────────
run_lm <- function(sub) {
  keep_cov <- COVARS[sapply(COVARS, function(v) {
    x <- sub[[v]]
    if (is.numeric(x)) sd(x, na.rm = TRUE) > 0 else length(unique(na.omit(x))) >= 2
  })]
  fmla <- as.formula(paste("proportion ~ BMI_group_num",
    if (length(keep_cov)) paste("+", paste(keep_cov, collapse = " + ")) else ""))
  m  <- lm(fmla, data = sub)
  ct <- summary(m)$coefficients
  data.frame(
    beta  = ct["BMI_group_num", "Estimate"],
    se    = ct["BMI_group_num", "Std. Error"],
    p_val = ct["BMI_group_num", "Pr(>|t|)"],
    n     = nrow(sub)
  )
}

lm_df <- prop_long %>%
  group_by(Depot, phenotype) %>%
  group_modify(~ run_lm(.x)) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p_val, method = "BH"))

lm_df$label <- with(lm_df, sprintf(
  "beta=%s\nadj.p=%s",
  formatC(beta, format = "f", digits = 4),
  formatC(p_adj, format = "e", digits = 1)
))

cat("\nLM results:\n")
print(as.data.frame(lm_df))

fwrite(lm_df, file.path(OUT_DIR, "09.SMC_phenotype_abundance_lm.tsv"),
       sep = "\t")
fwrite(prop_raw, file.path(OUT_DIR, "09.SMC_phenotype_proportions_per_donor.tsv"),
       sep = "\t")

# ─── Plot per depot: both phenotypes in one figure ───────────────────────────
# annotation: adj.p per phenotype, placed at top of each phenotype's boxplots
annot_df <- lm_df %>%
  mutate(label = paste0("beta=", formatC(beta, format = "f", digits = 3),
                        "\nadj.p=", formatC(p_adj, format = "e", digits = 2)))

# x-axis position per BMI group (numeric) for smooth — use BMI_group_num
# dodge width must match boxplot dodge
DODGE <- 0.75

for (dep in c("SAT", "VAT")) {
  sub   <- prop_long %>% filter(Depot == dep)
  lm_d  <- lm_df    %>% filter(Depot == dep)
  ann_d <- annot_df %>% filter(Depot == dep)

  # y position for annotation: per phenotype, above 95th pct of that phenotype's data
  ann_d <- ann_d %>%
    left_join(
      sub %>% group_by(phenotype) %>%
        summarise(ypos = quantile(proportion, 0.97, na.rm = TRUE) * 1.06,
                  .groups = "drop"),
      by = "phenotype"
    )

  p <- ggplot(sub, aes(x = BMI_group, y = proportion, fill = phenotype)) +
    geom_boxplot(
      position    = position_dodge(DODGE),
      outlier.size = 0.5, linewidth = 0.35, alpha = 0.8, width = 0.6
    ) +
    # LM trend lines per phenotype (use numeric x mapped to dodge positions)
    {
      lapply(seq_along(levels(sub$phenotype)), function(i) {
        ph    <- levels(sub$phenotype)[i]
        sig   <- lm_d$p_adj[lm_d$phenotype == ph] < 0.05
        lcol  <- PHENO_COLORS[[ph]]
        ltype <- if (sig) "solid" else "dashed"
        offset <- (i - 1.5) * DODGE / 2   # shift line to match dodge position
        geom_smooth(
          data        = sub %>% filter(phenotype == ph) %>%
                          mutate(x_dodge = BMI_group_num + offset),
          aes(x = x_dodge, y = proportion),
          method      = "lm", formula = y ~ x,
          color       = lcol, fill = lcol,
          linewidth   = 0.6, alpha = 0.12, linetype = ltype,
          inherit.aes = FALSE
        )
      })
    } +
    # adj.p annotation per phenotype (fixed x = 1 and 4 to avoid overlap)
    geom_text(
      data        = ann_d %>% mutate(xpos = ifelse(phenotype == "Contractile", 1, 4)),
      aes(x = xpos, y = ypos, label = label, color = phenotype),
      inherit.aes = FALSE,
      size = 2.3, hjust = 0.5, show.legend = FALSE
    ) +
    scale_fill_manual(values  = PHENO_COLORS, name = "SMC phenotype") +
    scale_color_manual(values = PHENO_COLORS) +
    scale_x_discrete(labels = function(x) gsub(" ", "\n", x)) +
    scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      expand = expansion(mult = c(0.02, 0.15))
    ) +
    labs(
      x     = NULL,
      y     = "Proportion of SMC cells",
      title = sprintf("%s — SMC phenotype proportion by BMI group", dep)
    ) +
    theme_classic(base_size = 9) +
    theme(
      plot.title      = element_text(size = 9, face = "bold", hjust = 0.5),
      axis.text.x     = element_text(size = 7),
      axis.text.y     = element_text(size = 7),
      axis.title.y    = element_text(size = 8),
      legend.position = "bottom",
      legend.title    = element_text(size = 7.5),
      legend.text     = element_text(size = 7),
      legend.key.size = unit(0.35, "cm")
    )

  out <- file.path(OUT_DIR, sprintf("09.SMC_phenotype_proportion_BMI_%s.pdf", dep))
  ggsave(out, p, width = 3.8, height = 3.4, device = "pdf")
  cat(sprintf("Saved: %s\n", basename(out)))
}

# ─── SAT vs VAT boxplot per phenotype (paired Wilcoxon) ──────────────────────
prop_wide <- prop_raw %>%
  select(donor_id, Depot, prop_contractile, prop_collagen) %>%
  filter(donor_id %in% (prop_raw %>%
    group_by(donor_id) %>% filter(n_distinct(Depot) == 2) %>% pull(donor_id) %>% unique()))

wilcox_stats <- lapply(c("Contractile" = "prop_contractile",
                          "Collagen/ECM" = "prop_collagen"), function(col) {
  sat <- prop_wide %>% filter(Depot == "SAT") %>% arrange(donor_id) %>% pull(col)
  vat <- prop_wide %>% filter(Depot == "VAT") %>% arrange(donor_id) %>% pull(col)
  wt  <- wilcox.test(sat, vat, paired = TRUE)
  data.frame(W = wt$statistic, p_value = wt$p.value,
             median_SAT = median(sat), median_VAT = median(vat))
}) %>% bind_rows(.id = "phenotype") %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))

cat("\nSAT vs VAT Wilcoxon results:\n")
print(as.data.frame(wilcox_stats))
fwrite(wilcox_stats, file.path(OUT_DIR, "09.SMC_phenotype_abundance_stats.tsv"), sep = "\t")

# Boxplot: SAT vs VAT per phenotype
prop_box <- prop_raw %>%
  select(donor_id, Depot, prop_contractile, prop_collagen) %>%
  tidyr::pivot_longer(cols = c(prop_contractile, prop_collagen),
                      names_to = "phenotype", values_to = "proportion") %>%
  mutate(phenotype = recode(phenotype,
    prop_contractile = "Contractile", prop_collagen = "Collagen/ECM"),
    phenotype = factor(phenotype, levels = c("Contractile", "Collagen/ECM")),
    Depot = factor(Depot, levels = c("SAT", "VAT")))

# annotation: exact p per phenotype (top of each facet)
annot_box <- wilcox_stats %>%
  mutate(
    label  = sprintf("p = %.2e", p_value),
    Depot  = "VAT",
    proportion = max(prop_box$proportion, na.rm = TRUE) * 1.06,
    phenotype = factor(phenotype, levels = c("Contractile", "Collagen/ECM"))
  )

p_box <- ggplot(prop_box, aes(x = Depot, y = proportion, fill = Depot)) +
  geom_boxplot(outlier.size = 0.5, linewidth = 0.35, alpha = 0.8, width = 0.55) +
  geom_text(data = annot_box,
            aes(x = 1.5, y = proportion, label = label),
            inherit.aes = FALSE, size = 2.5, color = "grey20", hjust = 0.5) +
  facet_wrap(~ phenotype, nrow = 1) +
  scale_fill_manual(values = DEPOT_COLORS) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0.02, 0.12))) +
  labs(x = NULL, y = "Proportion of SMC cells",
       title = "SMC phenotype proportion: SAT vs VAT") +
  theme_classic(base_size = 9) +
  theme(
    plot.title      = element_text(size = 9, face = "bold", hjust = 0.5),
    strip.text      = element_text(size = 8, face = "bold"),
    strip.background = element_blank(),
    legend.position = "none",
    axis.text       = element_text(size = 8)
  )

ggsave(file.path(OUT_DIR, "09.SMC_phenotype_abundance_boxplot.pdf"),
       p_box, width = 4, height = 3.5, device = "pdf")
cat("Saved: 09.SMC_phenotype_abundance_boxplot.pdf\n")

cat("\n[Done]\n")
