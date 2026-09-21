#!/usr/bin/env Rscript
# Purpose: Test donor-level associations between SP- and DP-ATM fractions in VAT.
# Input: SECTION4_METADATA.
# Outputs: Donor-level analysis data in RDS format and raw/adjusted model
# statistics in TSV format under SECTION4_OUTPUT.
# Method: Fractions are calculated relative to all VAT macrophages and tested
# with donor covariate adjustment after excluding Caucasian donors.
# Usage: Rscript sp_dp_fraction_association.R
library(readr)
library(dplyr)

# ── config ──────────────────────────────────────────────────
meta_path   <- Sys.getenv("SECTION4_METADATA")
output_dir <- Sys.getenv("SECTION4_OUTPUT")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

DEPOT_COL  <- "Depot"
DEPOT_VAL  <- "VAT"
DONOR_COL  <- "donor_id"
LV1_COL    <- "Annotation_level1"
LV2_COL    <- "Annotation_level2"

SP_LABEL   <- "TIM4+_ATM"
DP_LABEL   <- "TIM4+_CD11c+_ATM"
COV_COLS   <- c("Age", "Gender", "Race", "BMI", "Self_defined_DM")
MIN_N      <- 10

safe_as_numeric <- function(x) suppressWarnings(as.numeric(as.character(x)))

# ── load & filter VAT ────────────────────────────────────────
meda     <- read_tsv(meta_path, show_col_types = FALSE)
meda_vat <- meda %>% filter(.data[[DEPOT_COL]] == DEPOT_VAL)

# Donor covariates; Caucasian donors are excluded.
donor_covariates <- meda_vat %>%
  select(all_of(c(DONOR_COL, COV_COLS))) %>%
  group_by(.data[[DONOR_COL]]) %>%
  summarise(across(everything(), ~ first(na.omit(.x))), .groups = "drop") %>%
  filter(Race != "Caucasian") %>%
  mutate(
    Age             = safe_as_numeric(Age),
    BMI             = safe_as_numeric(BMI),
    Gender          = as.factor(Gender),
    Race            = as.factor(Race),
    Self_defined_DM = factor(
      if_else(Self_defined_DM == "T2DM", "T2DM", "Non_T2DM"),
      levels = c("Non_T2DM", "T2DM")
    )
  )

# ── cell counts ──────────────────────────────────────────────
macrophage_counts <- meda_vat %>%
  filter(.data[[LV1_COL]] == "Macrophage") %>%
  count(.data[[DONOR_COL]], name = "n_mac")

sp_counts <- meda_vat %>%
  filter(.data[[LV2_COL]] == SP_LABEL) %>%
  count(.data[[DONOR_COL]], name = "n_SP")

dp_counts <- meda_vat %>%
  filter(.data[[LV2_COL]] == DP_LABEL) %>%
  count(.data[[DONOR_COL]], name = "n_DP")

# Calculate SP- and DP-ATM fractions.
subtype_fractions <- macrophage_counts %>%
  left_join(sp_counts, by = DONOR_COL) %>%
  left_join(dp_counts, by = DONOR_COL) %>%
  mutate(across(c(n_SP, n_DP), ~ coalesce(.x, 0L))) %>%
  filter(n_SP >= MIN_N, n_DP >= MIN_N) %>%
  mutate(
    SP_frac_mac = n_SP / n_mac,
    DP_frac_mac = n_DP / n_mac
  )

# ── covariate join ───────────────────────────────────────────
analysis_data <- subtype_fractions %>%
  left_join(donor_covariates, by = DONOR_COL) %>%
  filter(complete.cases(.))

cat("[INFO] donors after filter:", nrow(analysis_data), "\n")

# Fit the raw model used for the displayed slope and the covariate-adjusted
# model used for inference.
raw_fit <- lm(SP_frac_mac ~ DP_frac_mac, data = analysis_data)
adjusted_fit <- lm(
  SP_frac_mac ~ DP_frac_mac + Age + Gender + Race + BMI + Self_defined_DM,
  data = analysis_data
)
raw_stat <- coef(summary(raw_fit))["DP_frac_mac", ]
adjusted_stat <- coef(summary(adjusted_fit))["DP_frac_mac", ]
model_statistics <- tibble(
  model = c("raw", "covariate_adjusted"),
  n = c(nobs(raw_fit), nobs(adjusted_fit)),
  beta = c(raw_stat[["Estimate"]], adjusted_stat[["Estimate"]]),
  standard_error = c(raw_stat[["Std. Error"]], adjusted_stat[["Std. Error"]]),
  p_value = c(raw_stat[["Pr(>|t|)"]], adjusted_stat[["Pr(>|t|)"]])
)

# Save analysis data and statistical results.
saveRDS(analysis_data, file.path(output_dir, "sp_dp_fraction_association.rds"))
write_tsv(
  model_statistics,
  file.path(output_dir, "sp_dp_fraction_association_statistics.tsv")
)
cat("[saved] sp_dp_fraction_association.rds and statistics TSV\n")
