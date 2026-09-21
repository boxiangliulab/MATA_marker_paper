#!/usr/bin/env Rscript
# Purpose: Test the association between donor DP-ATM status and T2D in VAT.
# Input: SECTION4_METADATA.
# Output: SECTION4_OUTPUT/dp_status_t2d_association.rds.
# Method: Donors with no DP ATMs are DP-absent, donors with at least 10 DP ATMs
# are DP-present, and intermediate donors are excluded. The unadjusted
# association is tested by chi-square; logistic regression adjusts for age,
# gender, race, and BMI.
# Usage: Rscript dp_status_t2d_association.R
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(broom)
})

# ── paths ─────────────────────────────────────────────────────
META    <- Sys.getenv("SECTION4_METADATA")
OUT_RDS <- file.path(Sys.getenv("SECTION4_OUTPUT"), "dp_status_t2d_association.rds")
dir.create(dirname(OUT_RDS), showWarnings = FALSE, recursive = TRUE)

# ── load & filter VAT ─────────────────────────────────────────
KEEP_COLS <- c("donor_id", "Depot", "Annotation_level2",
               "Age", "Gender", "Race", "BMI", "status",
               "Self_defined_DM", "b_FBG", "b_HbA1c",
               "b_Insulin (mU/L)", "Body Fat (%)")
meda <- fread(META, sep = "\t", data.table = FALSE, select = KEEP_COLS)
vat  <- meda %>% filter(Depot == "VAT")

# ── per-donor cell counts ─────────────────────────────────────
cnt_dp <- vat %>%
  filter(Annotation_level2 == "TIM4+_CD11c+_ATM") %>%
  count(donor_id, name = "n_DP")
cnt_sp <- vat %>%
  filter(Annotation_level2 == "TIM4+_ATM") %>%
  count(donor_id, name = "n_SP")

DONOR_META_COLS <- c(
  "donor_id", "Age", "Gender", "Race", "BMI", "status",
  "Self_defined_DM", "b_FBG", "b_HbA1c",
  "b_Insulin (mU/L)", "Body Fat (%)"
)
donor_meta <- vat %>%
  select(all_of(DONOR_META_COLS)) %>%
  distinct(donor_id, .keep_all = TRUE)

# ── assign DP/SP groups ───────────────────────────────────────
donor_df <- donor_meta %>%
  left_join(cnt_dp, by = "donor_id") %>%
  left_join(cnt_sp, by = "donor_id") %>%
  mutate(
    n_DP = replace_na(n_DP, 0),
    n_SP = replace_na(n_SP, 0),
    group_DP = case_when(
      n_DP >= 10 ~ "DP present",
      n_DP == 0  ~ "DP absent",
      TRUE       ~ NA_character_
    ),
    group_SP = case_when(
      n_SP >= 10 ~ "SP present",
      n_SP == 0  ~ "SP absent",
      TRUE       ~ NA_character_
    ),
    Self_defined_DM = factor(
      ifelse(Self_defined_DM == "T2DM", "T2DM", "Non_T2DM"),
      levels = c("Non_T2DM", "T2DM")
    ),
    DM_binary = as.character(Self_defined_DM)
  ) %>%
  filter(Race != "Caucasian")

cat(sprintf("total donors: %d\n", nrow(donor_df)))
print(table(donor_df$group_DP, useNA = "ifany"))
print(table(donor_df$group_SP, useNA = "ifany"))
print(table(donor_df$DM_binary, useNA = "ifany"))

# ── helper: Chi-square ────────────────────────────────────────
run_test <- function(dm_vec, group_vec) {
  df <- data.frame(DM = dm_vec, Group = group_vec) %>%
    filter(!is.na(DM), !is.na(Group))
  tbl   <- table(DM = df$DM, Group = df$Group)
  p_val <- suppressWarnings(chisq.test(tbl)$p.value)
  list(
    result = data.frame(
      test    = "Chi-square",
      p_value = p_val,
      sig     = case_when(
        p_val < 0.001 ~ "***",
        p_val < 0.01  ~ "**",
        p_val < 0.05  ~ "*",
        TRUE          ~ "ns"
      ),
      n_total = sum(tbl),
      stringsAsFactors = FALSE
    ),
    contingency = tbl
  )
}

# ── Chi-square: DP ────────────────────────────────────────────
dp_test_df <- donor_df %>% filter(!is.na(group_DP))
dp_res     <- run_test(dp_test_df$DM_binary, dp_test_df$group_DP)
cat("\n[DP] Contingency table:\n"); print(dp_res$contingency)
cat("\n[DP] Test result:\n");       print(dp_res$result)

# ── Chi-square: SP ────────────────────────────────────────────
sp_test_df <- donor_df %>% filter(!is.na(group_SP))
sp_res     <- run_test(sp_test_df$DM_binary, sp_test_df$group_SP)
cat("\n[SP] Contingency table:\n"); print(sp_res$contingency)
cat("\n[SP] Test result:\n");       print(sp_res$result)

# Covariate-adjusted DP-status association reported in the manuscript.
dp_model_df <- donor_df %>%
  filter(!is.na(group_DP), !is.na(DM_binary)) %>%
  mutate(
    group_DP = factor(group_DP, levels = c("DP absent", "DP present")),
    T2D = as.integer(DM_binary == "T2DM"),
    Age = suppressWarnings(as.numeric(as.character(Age))),
    BMI = suppressWarnings(as.numeric(as.character(BMI))),
    Gender = factor(Gender),
    Race = factor(Race)
  ) %>%
  drop_na(T2D, group_DP, Age, Gender, Race, BMI)

dp_adjusted_fit <- glm(
  T2D ~ group_DP + Age + Gender + Race + BMI,
  data = dp_model_df,
  family = binomial
)
dp_adjusted_result <- tidy(dp_adjusted_fit, conf.int = TRUE) %>%
  filter(term == "group_DPDP present") %>%
  transmute(
    model = "age_gender_race_bmi_adjusted",
    n = nrow(dp_model_df),
    odds_ratio = exp(estimate),
    conf_low = exp(conf.low),
    conf_high = exp(conf.high),
    p_value = p.value
  )
cat("\n[DP] Covariate-adjusted logistic regression:\n")
print(dp_adjusted_result)

# ── save RDS ──────────────────────────────────────────────────
obj <- list(
  generated_at   = as.character(Sys.time()),
  meta_path      = META,
  donor_df       = donor_df,
  chisq_res_dp   = dp_res$result,
  contingency_dp = dp_res$contingency,
  chisq_res_sp   = sp_res$result,
  contingency_sp = sp_res$contingency,
  adjusted_logistic_dp = dp_adjusted_result
)
saveRDS(obj, OUT_RDS)
fwrite(
  as.data.table(dp_adjusted_result),
  file.path(Sys.getenv("SECTION4_OUTPUT"), "dp_status_t2d_adjusted_logistic.tsv"),
  sep = "\t"
)
cat(sprintf("\n[SAVED] %s\n", OUT_RDS))
