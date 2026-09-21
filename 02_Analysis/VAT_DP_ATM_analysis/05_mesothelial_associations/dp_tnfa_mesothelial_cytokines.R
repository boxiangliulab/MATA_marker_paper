#!/usr/bin/env Rscript
# Purpose: Test associations between the DP-ATM TNFA/NFKB module score and
# mesothelial IL6 or CXCL8 expression in VAT.
# Inputs: SECTION4_DP_MODULE_SCORES, SECTION4_MESOTHELIAL_PSEUDOBULK, and
# SECTION4_METADATA.
# Output: SECTION4_OUTPUT/dp_tnfa_mesothelial_cytokines.rds.
# Method: Partial correlations and covariate-adjusted linear models among donors
# with at least 10 DP ATMs and 10 mesothelial cells.
# Usage: Rscript dp_tnfa_mesothelial_cytokines.R
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
})

# ── paths ──────────────────────────────────────────────────────
dp_ms_path   <- Sys.getenv("SECTION4_DP_MODULE_SCORES")
meso_pb_path <- Sys.getenv("SECTION4_MESOTHELIAL_PSEUDOBULK")
meta_path    <- Sys.getenv("SECTION4_METADATA")
OUT_RDS      <- file.path(Sys.getenv("SECTION4_OUTPUT"), "dp_tnfa_mesothelial_cytokines.rds")
dir.create(dirname(OUT_RDS), showWarnings = FALSE, recursive = TRUE)

DEPOT_COL    <- "Depot"
DEPOT_VAL    <- "VAT"
DONOR_COL    <- "donor_id"
MS_DONOR_COL <- "donor"
LV2_COL      <- "Annotation_level2"
MESO_LABEL   <- "Mesothelial"
DP_LABEL     <- "TIM4+_CD11c+_ATM"
DP_MS_COL    <- "MS_TNFA_NFKB"
TARGET_GENES <- c("IL6", "CXCL8")
COV_COLS     <- c("Age", "Gender", "Race", "BMI", "Self_defined_DM")
MIN_CELLS    <- 10
MIN_N        <- 10

safe_as_numeric <- function(x) suppressWarnings(as.numeric(as.character(x)))

# ── metadata -> covariates + cell filter ───────────────────────
meda <- read_tsv(meta_path, show_col_types = FALSE)

donors_dp_pass <- meda %>%
  filter(.data[[DEPOT_COL]] == DEPOT_VAL, .data[[LV2_COL]] == DP_LABEL) %>%
  count(.data[[DONOR_COL]], name = "n_dp_cells") %>%
  filter(n_dp_cells >= MIN_CELLS) %>%
  pull(.data[[DONOR_COL]])

donors_meso_pass <- meda %>%
  filter(.data[[DEPOT_COL]] == DEPOT_VAL, .data[[LV2_COL]] == MESO_LABEL) %>%
  count(.data[[DONOR_COL]], name = "n_meso_cells") %>%
  filter(n_meso_cells >= MIN_CELLS) %>%
  pull(.data[[DONOR_COL]])

donors_cells_pass <- intersect(donors_dp_pass, donors_meso_pass)

cat("[INFO] Donors >= ", MIN_CELLS, " DP cells: ",   length(donors_dp_pass),    "\n", sep = "")
cat("[INFO] Donors >= ", MIN_CELLS, " Meso cells: ", length(donors_meso_pass),  "\n", sep = "")
cat("[INFO] Both cell filters: ",                    length(donors_cells_pass), "\n", sep = "")

cov_df <- meda %>%
  filter(.data[[DEPOT_COL]] == DEPOT_VAL) %>%
  select(all_of(c(DONOR_COL, COV_COLS))) %>%
  group_by(.data[[DONOR_COL]]) %>%
  summarise(across(everything(), ~ first(na.omit(.x))), .groups = "drop") %>%
  filter(.data[[DONOR_COL]] %in% donors_cells_pass) %>%
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

# ── DP module score ────────────────────────────────────────────
dp_df <- read_tsv(dp_ms_path, show_col_types = FALSE) %>%
  select(all_of(c(MS_DONOR_COL, DP_MS_COL))) %>%
  rename(!!DONOR_COL := all_of(MS_DONOR_COL)) %>%
  distinct(.data[[DONOR_COL]], .keep_all = TRUE)

# ── Mesothelial pseudobulk ────────────────────────────────────
meso_raw  <- read_tsv(meso_pb_path, show_col_types = FALSE)
gene_col  <- if ("Gene" %in% colnames(meso_raw)) "Gene" else colnames(meso_raw)[1]
meso_mat  <- meso_raw %>%
  mutate(across(-all_of(gene_col), as.numeric)) %>%
  column_to_rownames(gene_col) %>%
  as.matrix()

TARGET_GENES <- intersect(TARGET_GENES, rownames(meso_mat))

meso_df <- as.data.frame(t(meso_mat[TARGET_GENES, , drop = FALSE])) %>%
  rownames_to_column(DONOR_COL)

# ── align donors ──────────────────────────────────────────────
donors_common <- Reduce(intersect, list(
  dp_df[[DONOR_COL]],
  meso_df[[DONOR_COL]],
  cov_df[[DONOR_COL]]
)) %>% sort()

cat("[INFO] Common donors after all filters:", length(donors_common), "\n")

dp_sub   <- dp_df   %>% filter(.data[[DONOR_COL]] %in% donors_common) %>%
              arrange(match(.data[[DONOR_COL]], donors_common))
meso_sub <- meso_df %>% filter(.data[[DONOR_COL]] %in% donors_common) %>%
              arrange(match(.data[[DONOR_COL]], donors_common))
cov_sub  <- cov_df  %>% filter(.data[[DONOR_COL]] %in% donors_common) %>%
              arrange(match(.data[[DONOR_COL]], donors_common)) %>%
              select(all_of(COV_COLS))

# ── partial correlation ────────────────────────────────────────
pcor_residual <- function(x, y, cov_df) {
  df <- data.frame(x = as.numeric(x), y = as.numeric(y), cov_df)
  df <- df[complete.cases(df), , drop = FALSE]
  n  <- nrow(df)
  if (n < MIN_N) return(list(r = NA_real_, p = NA_real_, n = n, beta = NA_real_, df_resid = NA_integer_))
  X        <- stats::model.matrix(~ ., data = df[, -(1:2), drop = FALSE])
  e_x      <- resid(lm(df$x ~ X - 1))
  e_y      <- resid(lm(df$y ~ X - 1))
  df_resid <- n - ncol(X)
  if (df_resid < 2) return(list(r = NA_real_, p = NA_real_, n = n, beta = NA_real_, df_resid = df_resid))
  r        <- unname(cor(e_x, e_y))
  tstat    <- r * sqrt(df_resid / max(1e-12, 1 - r^2))
  p        <- 2 * stats::pt(-abs(tstat), df = df_resid)
  beta_raw <- unname(coef(lm(df$y ~ df$x))["df$x"])
  list(r = r, p = p, n = n, beta = beta_raw, df_resid = df_resid)
}

run_pcor <- function(remove_zero) {
  lapply(TARGET_GENES, function(g) {
    keep <- if (remove_zero) !is.na(meso_sub[[g]]) & meso_sub[[g]] != 0
            else             !is.na(meso_sub[[g]])
    tmp <- pcor_residual(dp_sub[[DP_MS_COL]][keep], meso_sub[[g]][keep], cov_sub[keep, ])
    data.frame(
      dp_score     = DP_MS_COL,
      meso_gene    = g,
      r            = tmp$r,
      beta         = tmp$beta,
      p            = tmp$p,
      n            = tmp$n,
      df_resid     = tmp$df_resid,
      zero_removed = remove_zero,
      stringsAsFactors = FALSE
    )
  }) %>% bind_rows()
}

res_nozero   <- run_pcor(remove_zero = TRUE)
res_withzero <- run_pcor(remove_zero = FALSE)

cat("--- zero removed ---\n"); print(res_nozero)
cat("--- zero kept ---\n");    print(res_withzero)

# ── save ──────────────────────────────────────────────────────
obj <- list(
  generated_at  = as.character(Sys.time()),
  dp_ms_col     = DP_MS_COL,
  target_genes  = TARGET_GENES,
  donors_common = donors_common,
  dp_sub        = dp_sub,
  meso_sub      = meso_sub,
  cov_sub       = cov_sub,
  res_nozero    = res_nozero,
  res_withzero  = res_withzero
)

saveRDS(obj, OUT_RDS)
cat(sprintf("\n[SAVED] %s\n", OUT_RDS))
