library(tidyverse)
library(ggplot2)
library(cowplot)

# Project root is 7 levels up from this script
args0 <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args0[grep("--file=", args0)])
project_root <- normalizePath(file.path(dirname(normalizePath(script_path)),
                                        paste(rep("..", 7), collapse = "/")))

OUT_DIR <- file.path(project_root, "data/30.adipogenesis")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---------- 0. Data ----------------------------------------------------------
meta <- read.table(
  file.path(project_root, "data/08.integration_snRNA_object/01.freeze_version/metadata_after_libraryQC_after_ctQC_freeze.tsv"),
  sep = "\t", header = TRUE
)

# ---------- 1. Subset Asian donors & target cell types -----------------------
CELL_TYPES  <- c("Early_preadipocyte", "Committed_preadipocyte", "Adipocyte")
ASIAN_RACES <- c("Chinese", "Malay", "Indian")

meta_asian <- meta %>%
  filter(Race %in% ASIAN_RACES) %>%
  filter(Annotation_level2 %in% CELL_TYPES)

# ---------- 2. Calculate cell proportions ------------------------------------
# Total nuclei per donor-depot (all cell types, before subsetting)
total_nuclei <- meta %>%
  filter(Race %in% ASIAN_RACES) %>%
  group_by(donor_id, Depot) %>%
  summarise(total = n(), .groups = "drop")

# Nuclei per donor-depot-celltype
ct_counts <- meta_asian %>%
  group_by(donor_id, Depot, Annotation_level2) %>%
  summarise(ct_n = n(), .groups = "drop")

prop_df <- ct_counts %>%
  left_join(total_nuclei, by = c("donor_id", "Depot")) %>%
  # QC: exclude donor-depot with < 100 total nuclei
  filter(total >= 100) %>%
  # QC: exclude cell subtypes with < 10 nuclei in that donor-depot
  filter(ct_n >= 10) %>%
  mutate(proportion       = ct_n / total,
         log10_proportion = log10(proportion))

# ---------- 3. Merge donor-level covariates ----------------------------------
donor_meta <- meta_asian %>%
  select(donor_id, BMI, Age, Gender, Diabetic, Race) %>%
  distinct(donor_id, .keep_all = TRUE) %>%
  mutate(T2DM = ifelse(Diabetic == "T2DM", 1, 0))

prop_df <- prop_df %>%
  left_join(donor_meta, by = "donor_id")

# ---------- 4. Asian BMI groups ----------------------------------------------
prop_df <- prop_df %>%
  mutate(BMI_group = case_when(
    BMI <  27.5 ~ "Non-obese",
    BMI <  32.5 ~ "Obese",
    BMI <  37.5 ~ "Severely obese",
    TRUE         ~ "Morbidly obese"
  )) %>%
  mutate(
    BMI_group = factor(BMI_group, levels = c(
      "Non-obese", "Obese", "Severely obese", "Morbidly obese"
    )),
    BMI_group_num = as.numeric(BMI_group)
  )

# ---------- 5. Plot ----------------------------------------------------------
GROUP_COLORS <- c(
  "Non-obese"      = "#efedf5",
  "Obese"          = "#bcbddc",
  "Severely obese" = "#807dba",
  "Morbidly obese" = "#54278f"
)

CELL_LABELS <- c(
  "Early_preadipocyte"     = "Early preadipocyte",
  "Committed_preadipocyte" = "Committed preadipocyte",
  "Adipocyte"              = "Adipocyte"
)

prop_df <- prop_df %>%
  mutate(Annotation_level2 = factor(
    Annotation_level2,
    levels = c("Early_preadipocyte", "Committed_preadipocyte", "Adipocyte")
  ))

# Multivariable linear regression: log10(proportion) ~ BMI_group_num + Age + Gender + T2DM + Race
# Covariates with < 2 levels or zero SD within the stratum are dropped automatically
fit_lm_safe <- function(df) {
  covars <- c("Age", "Gender", "T2DM", "Race")
  keep <- covars[sapply(covars, function(v) {
    x <- df[[v]]
    if (is.numeric(x)) sd(x, na.rm = TRUE) > 0 else length(unique(na.omit(x))) >= 2
  })]
  fmla <- as.formula(paste("log10_proportion ~ BMI_group_num +", paste(keep, collapse = " + ")))
  lm(fmla, data = df)
}

fmt_sci <- function(p) formatC(p, format = "e", digits = 2)

stat_labels <- prop_df %>%
  group_by(Depot, Annotation_level2) %>%
  group_modify(~ {
    # LM: BMI_group_num + covariates
    m        <- fit_lm_safe(.x)
    coef_tbl <- summary(m)$coefficients
    lm_beta  <- coef_tbl["BMI_group_num", "Estimate"]
    lm_p     <- coef_tbl["BMI_group_num", "Pr(>|t|)"]
    # Spearman: BMI_group_num vs log10_proportion
    sp       <- cor.test(.x$BMI_group_num, .x$log10_proportion, method = "spearman", exact = FALSE)
    data.frame(
      lm_beta = lm_beta,
      lm_p    = lm_p,
      sp_rho  = sp$estimate,
      sp_p    = sp$p.value
    )
  }) %>%
  ungroup() %>%
  # BH correction across all 6 tests (3 cell types × 2 depots)
  mutate(
    lm_p_adj = p.adjust(lm_p, method = "BH"),
    sp_p_adj = p.adjust(sp_p, method = "BH"),
    label = paste0(
      "beta = ", round(lm_beta, 4), ", p.adj = ", fmt_sci(lm_p_adj), " (lm)\n",
      "Spearman R = ", round(sp_rho, 3), ", p.adj = ", fmt_sci(sp_p_adj)
    )
  )

cat("\n=== LM results with BH correction (6 tests) ===\n")
print(stat_labels %>% select(Depot, Annotation_level2, lm_beta, lm_p, lm_p_adj, sp_rho, sp_p, sp_p_adj))

# N per BMI group per depot x celltype
n_labels <- prop_df %>%
  group_by(Depot, Annotation_level2, BMI_group) %>%
  summarise(n = n(), .groups = "drop")

# significance flag based on BH-adjusted p
sig_dep_ct <- stat_labels %>%
  select(Depot, Annotation_level2, lm_p_adj) %>%
  mutate(trend = ifelse(lm_p_adj < 0.05, "sig", "ns"))
prop_df <- prop_df %>% left_join(sig_dep_ct %>% select(Depot, Annotation_level2, trend),
                                  by = c("Depot", "Annotation_level2"))

plot_one <- function(depot_name) {
  dat   <- prop_df    %>% filter(Depot == depot_name)
  stats <- stat_labels %>% filter(Depot == depot_name)
  n_lab <- n_labels   %>% filter(Depot == depot_name)

  ggplot(dat, aes(x = BMI_group, y = log10_proportion, fill = BMI_group)) +
    geom_boxplot(outlier.size = 0.5, width = 0.35, color = "grey30") +
    geom_jitter(width = 0.1, size = 0.3, color = "grey40", alpha = 0.5) +
    geom_smooth(
      data = dat %>% filter(trend == "sig"),
      aes(x = BMI_group_num, y = log10_proportion),
      method = "lm", formula = y ~ x,
      color = "red", fill = "red", alpha = 0.15,
      linewidth = 0.5, inherit.aes = FALSE
    ) +
    geom_smooth(
      data = dat %>% filter(trend == "ns"),
      aes(x = BMI_group_num, y = log10_proportion),
      method = "lm", formula = y ~ x,
      color = "grey50", fill = "grey70", linetype = "dashed",
      linewidth = 0.5, alpha = 0.15, inherit.aes = FALSE
    ) +
    geom_text(
      data = n_lab,
      aes(x = BMI_group, y = -Inf, label = paste0("n=", n)),
      inherit.aes = FALSE,
      vjust = -0.3, size = 2.2, color = "grey30"
    ) +
    facet_wrap(~ Annotation_level2, nrow = 1,
               labeller = labeller(Annotation_level2 = CELL_LABELS)) +
    geom_text(
      data = stats,
      aes(x = 2.5, y = Inf, label = label),
      inherit.aes = FALSE,
      vjust = 1.2, hjust = 0.5, size = 2.3
    ) +
    scale_fill_manual(values = GROUP_COLORS, name = "BMI group") +
    labs(x = NULL, y = expression(log[10](proportion)), title = depot_name) +
    theme_classic(base_size = 9, base_family = "Helvetica") +
    theme(
      legend.position  = "none",
      axis.text.x      = element_text(size = 7, angle = 90, vjust = 0.5, hjust = 1),
      strip.background = element_blank(),
      strip.text       = element_text(size = 8, face = "bold"),
      plot.title       = element_text(size = 9, face = "bold"),
      aspect.ratio     = 1.5
    )
}

p_sat <- plot_one("SAT")
p_vat <- plot_one("VAT")

# Extract legend from a separate copy with legend enabled
legend_plot <- p_sat +
  theme(
    legend.position  = "right",
    legend.key.size  = unit(0.25, "cm"),
    legend.text      = element_text(size = 6),
    legend.title     = element_text(size = 6, face = "bold")
  )
shared_legend <- get_legend(legend_plot)

# Panels without legend (6 wide) + legend strip on right (~1.2 wide)
p_panels   <- plot_grid(p_sat, p_vat, nrow = 2)
p_combined <- plot_grid(p_panels, shared_legend, nrow = 1, rel_widths = c(6, 1.2))

out_file <- file.path(OUT_DIR, "01.cell_proportion_BMI_asian.pdf")
ggsave(filename = out_file, plot = p_combined, width = 5.5, height = 6)
message("Done: ", out_file)
