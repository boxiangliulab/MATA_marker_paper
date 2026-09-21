#!/usr/bin/env Rscript
# Purpose: Test mesothelial fibrosis-gene expression by donor DP-ATM status.
# Inputs: SECTION4_MESOTHELIAL_PSEUDOBULK and SECTION4_VAT_DONOR_GROUPS.
# Output: SECTION4_OUTPUT/mesothelial_fibrosis_by_dp_status.rds.
# Method: Linear models adjusted for age, gender, race, BMI, and diabetes status.
# Usage: Rscript mesothelial_fibrosis_by_dp_status.R

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(broom)
})

# ── Paths ──────────────────────────────────────────────────────────────────────
EXPR_FILE <- Sys.getenv("SECTION4_MESOTHELIAL_PSEUDOBULK")
GRP_FILE  <- Sys.getenv("SECTION4_VAT_DONOR_GROUPS")
OUT_RDS   <- file.path(Sys.getenv("SECTION4_OUTPUT"), "mesothelial_fibrosis_by_dp_status.rds")
dir.create(dirname(OUT_RDS), showWarnings = FALSE, recursive = TRUE)

GENES   <- c("COL1A1", "CCN2", "SERPINE1", "TGFBI", "FN1", "TNC")
DP_LEVS <- c("DP_absent", "DP_present")

sig_label <- function(p) {
  case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "ns")
}

# ── Load & merge ───────────────────────────────────────────────────────────────
cat("Loading pseudobulk expression...\n")
expr_raw <- read.table(EXPR_FILE, sep = "\t", header = TRUE,
                       row.names = 1, check.names = FALSE)

expr_df <- as.data.frame(t(expr_raw[GENES, , drop = FALSE])) %>%
  rownames_to_column("donor_id")

cat("Loading donor group info...\n")
grp <- read.table(GRP_FILE, sep = "\t", header = TRUE,
                  check.names = FALSE, na.strings = c("NA", "", "x"))

df <- expr_df %>%
  inner_join(grp, by = "donor_id") %>%
  mutate(
    group_DP        = factor(group_DP, levels = DP_LEVS),
    Gender          = factor(Gender),
    Race            = factor(Race),
    Self_defined_DM = factor(
      ifelse(Self_defined_DM == "T2DM", "T2DM", "Non_T2DM"),
      levels = c("Non_T2DM", "T2DM")
    ),
    BMI             = as.numeric(BMI),
    Age             = as.numeric(Age)
  ) %>%
  filter(Race != "Caucasian")

df_dp <- df %>% filter(!is.na(group_DP))

cat(sprintf("Working dataset: n=%d  (DP_present=%d, DP_absent=%d)\n",
    nrow(df_dp),
    sum(df_dp$group_DP == "DP_present"),
    sum(df_dp$group_DP == "DP_absent")))

# ── Linear regression (M2 only) ───────────────────────────────────────────────
cat("Running linear regression...\n")

run_lm <- function(data, gene, rhs, model_name) {
  dat <- data %>%
    select(all_of(c(gene, "group_DP", "Age", "Gender", "Race",
                    "BMI", "Self_defined_DM"))) %>%
    drop_na()
  fit <- lm(as.formula(paste(gene, "~", rhs)), data = dat)
  tidy(fit, conf.int = TRUE) %>%
    filter(term == "group_DPDP_present") %>%
    mutate(gene = gene, model = model_name, n = nrow(dat),
           sig = sig_label(p.value), .before = 1) %>%
    select(gene, model, n, estimate, std.error,
           conf.low, conf.high, statistic, p.value, sig)
}

lm_results <- list()
for (g in GENES) {
  res <- run_lm(df_dp, g, "group_DP + Age + Gender + Race + BMI + Self_defined_DM", "M2")
  lm_results[[g]] <- res
}

lm_df <- bind_rows(lm_results) %>%
  mutate(
    p.adj     = p.adjust(p.value, method = "BH"),
    sig       = sig_label(p.adj)
  )

# ── annotation lookup ─────────────────────────────────────────────────────────
ann_lookup <- lm_df %>%
  select(gene, model, estimate, p.value, p.adj, sig) %>%
  pivot_wider(names_from = model,
              values_from = c(estimate, p.value, p.adj, sig),
              names_sep   = "_") %>%
  mutate(sig_display = sig_M2)

# ── Save RDS ───────────────────────────────────────────────────────────────────
obj <- list(
  generated_at = as.character(Sys.time()),
  expr_file    = EXPR_FILE,
  grp_file     = GRP_FILE,
  genes        = GENES,
  dp_levs      = DP_LEVS,
  df_dp        = df_dp,
  lm_df        = lm_df,
  ann_lookup   = ann_lookup
)

saveRDS(obj, OUT_RDS)
cat(sprintf("\n[SAVED] %s\n", OUT_RDS))
