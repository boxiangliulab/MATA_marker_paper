#!/usr/bin/env Rscript
# Purpose: Test associations of DP ATM abundance and DP presence with baseline
# and three-year HbA1c among non-Caucasian VAT donors.
# Inputs: SECTION4_METADATA.
# Output: SECTION4_OUTPUT/dp_status_hba1c_analysis.tsv.
# Models: HbA1c outcomes are tested against either log10 DP ATM proportion or
# binary DP status with the covariates specified below.
# Usage: Rscript dp_status_hba1c_analysis.R
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

META <- Sys.getenv("SECTION4_METADATA")
OUT <- file.path(Sys.getenv("SECTION4_OUTPUT"), "dp_status_hba1c_analysis.tsv")
dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)
three_year_col <- "X3yr_HbA1c"
surgery_col <- "Surg"
selected_cols <- unique(c(
  "donor_id", "Depot", "Annotation_level1", "Annotation_level2", "b_HbA1c",
  "Age", "Gender", "BMI", "Race", three_year_col, surgery_col
))

m <- fread(META, select = selected_cols, data.table = FALSE)
setnames(m, three_year_col, "H3")
setnames(m, surgery_col, "Surg")

# The frozen metadata correction for donor H734 is retained from the analysis.
m <- m %>% mutate(b_HbA1c = ifelse(donor_id == "H734", "5.8", b_HbA1c))
vat <- m %>% filter(Depot == "VAT")
mac <- vat %>% filter(Annotation_level1 == "Macrophage") %>% count(donor_id, name = "n_mac")
dp <- vat %>% filter(Annotation_level2 == "TIM4+_CD11c+_ATM") %>% count(donor_id, name = "n_dp")

meta <- vat %>%
  select(donor_id, b_HbA1c, H3, Age, Gender, BMI, Race, Surg) %>%
  distinct(donor_id, .keep_all = TRUE) %>%
  mutate(
    across(c(b_HbA1c, H3, Age, BMI), ~ suppressWarnings(as.numeric(na_if(as.character(.), "x")))),
    Gender = factor(Gender, levels = c("M", "F")),
    Race = factor(Race),
    Surg_group = factor(case_when(
      grepl("Sleeve", Surg, TRUE) ~ "Sleeve",
      grepl("RYGB", Surg, TRUE) ~ "RYGB",
      grepl("OAGB", Surg, TRUE) ~ "OAGB",
      TRUE ~ "Other"
    ))
  )

d <- mac %>%
  left_join(dp, by = "donor_id") %>%
  mutate(
    n_dp = tidyr::replace_na(n_dp, 0L),
    dp_prop = n_dp / n_mac,
    dp_bin = case_when(n_dp >= 10 ~ 1L, n_dp == 0 ~ 0L, TRUE ~ NA_integer_)
  ) %>%
  left_join(meta, by = "donor_id") %>%
  filter(Race != "Caucasian")

fitstat <- function(df, pred, out, cov) {
  df <- df %>% filter(!is.na(.data[[out]]), !is.na(Age), !is.na(Gender), !is.na(BMI), !is.na(Race))
  if (grepl("Surg|b_HbA1c", paste(cov, collapse = ""))) {
    df <- df %>% filter(!is.na(b_HbA1c), !is.na(Surg_group))
  }
  raw_formula <- as.formula(paste(out, "~", pred))
  adjusted_formula <- as.formula(paste(out, "~", pred, "+", paste(cov, collapse = "+")))
  raw_model <- lm(raw_formula, df)
  list(
    n = nrow(df),
    beta = unname(coef(raw_model)[pred]),
    p_raw = unname(coef(summary(raw_model))[pred, 4]),
    p = unname(coef(summary(lm(adjusted_formula, df)))[pred, 4])
  )
}
COVb <- c("Age", "Gender", "BMI", "Race")
COV3 <- c("b_HbA1c", "Age", "Gender", "BMI", "Race", "Surg_group")
pA <- d %>% filter(n_dp >= 10) %>% mutate(log_dp = log10(dp_prop))
pB <- d %>% filter(!is.na(dp_bin))
a1 <- fitstat(pA, "log_dp", "b_HbA1c", COVb)
a2 <- fitstat(pB, "dp_bin", "b_HbA1c", COVb)
b1 <- fitstat(pA, "log_dp", "H3", COV3)
b2 <- fitstat(pB, "dp_bin", "H3", COV3)

results <- bind_rows(
  data.frame(outcome = "baseline_HbA1c", exposure = "log10_DP_fraction", as.data.frame(a1)),
  data.frame(outcome = "baseline_HbA1c", exposure = "DP_presence", as.data.frame(a2)),
  data.frame(outcome = "three_year_HbA1c", exposure = "log10_DP_fraction", as.data.frame(b1)),
  data.frame(outcome = "three_year_HbA1c", exposure = "DP_presence", as.data.frame(b2))
)
dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)
fwrite(results, OUT, sep = "\t")
message("Saved: ", OUT)
