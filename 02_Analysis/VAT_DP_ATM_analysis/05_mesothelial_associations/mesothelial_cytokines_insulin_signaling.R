#!/usr/bin/env Rscript
# Purpose: Test associations between mesothelial IL6 or CXCL8 expression and
# mesothelial insulin-signalling module scores in VAT.
# Inputs: SECTION4_MESOTHELIAL_PSEUDOBULK,
# SECTION4_MESOTHELIAL_MODULE_SCORES, and SECTION4_METADATA.
# Output: SECTION4_OUTPUT/mesothelial_cytokines_insulin_signaling.rds.
# Method: Partial correlations and covariate-adjusted linear models among donors
# passing mesothelial and DP-ATM cell-count thresholds.
# Usage: Rscript mesothelial_cytokines_insulin_signaling.R
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
})

# ── paths ──────────────────────────────────────────────────────
meso_pb_path <- Sys.getenv("SECTION4_MESOTHELIAL_PSEUDOBULK")
meso_ms_path <- Sys.getenv("SECTION4_MESOTHELIAL_MODULE_SCORES")
meta_path    <- Sys.getenv("SECTION4_METADATA")

OUT_RDS <- file.path(Sys.getenv("SECTION4_OUTPUT"), "mesothelial_cytokines_insulin_signaling.rds")

DEPOT_COL    <- "Depot"
DEPOT_VAL    <- "VAT"
DONOR_COL    <- "donor_id"
MS_DONOR_COL <- "donor"
MESO_LABEL   <- "Mesothelial"
DP_LABEL     <- "TIM4+_CD11c+_ATM"
LV2_COL      <- "Annotation_level2"
TARGET_GENES <- c("IL6", "CXCL8")
MESO_MS_COL  <- "MS_REACTOME_IR_CASCADE"
COV_COLS     <- c("Age", "Gender", "Race", "BMI", "Self_defined_DM")
MIN_CELLS    <- 10
MIN_N        <- 10

safe_as_numeric <- function(x) suppressWarnings(as.numeric(as.character(x)))

# ── metadata -> covariates + Meso/DP cell filters ──────────────
meda <- read_tsv(meta_path, show_col_types = FALSE)

donors_meso_pass <- meda %>%
  filter(.data[[DEPOT_COL]] == DEPOT_VAL, .data[[LV2_COL]] == MESO_LABEL) %>%
  count(.data[[DONOR_COL]], name = "n_meso_cells") %>%
  filter(n_meso_cells >= MIN_CELLS) %>%
  pull(.data[[DONOR_COL]])

donors_dp_pass <- meda %>%
  filter(.data[[DEPOT_COL]] == DEPOT_VAL, .data[[LV2_COL]] == DP_LABEL) %>%
  count(.data[[DONOR_COL]], name = "n_dp_cells") %>%
  filter(n_dp_cells >= MIN_CELLS) %>%
  pull(.data[[DONOR_COL]])

donors_cell_pass <- intersect(donors_meso_pass, donors_dp_pass) %>% sort()

cat("[INFO] Donors >=", MIN_CELLS, "Meso cells:", length(donors_meso_pass), "\n")
cat("[INFO] Donors >=", MIN_CELLS, "DP cells:", length(donors_dp_pass), "\n")
cat("[INFO] Donors >=", MIN_CELLS, "Meso cells AND >=", MIN_CELLS, "DP cells:", length(donors_cell_pass), "\n")

cov_df <- meda %>%
  filter(.data[[DEPOT_COL]] == DEPOT_VAL) %>%
  dplyr::select(all_of(c(DONOR_COL, COV_COLS))) %>%
  group_by(.data[[DONOR_COL]]) %>%
  summarise(across(everything(), ~ first(na.omit(.x))), .groups = "drop") %>%
  filter(.data[[DONOR_COL]] %in% donors_cell_pass) %>%
  filter(Race != "Caucasian") %>%
  mutate(
    Age             = safe_as_numeric(Age),
    BMI             = safe_as_numeric(BMI),
    Gender          = as.factor(Gender),
    Race            = as.factor(Race),
    Self_defined_DM = factor(if_else(Self_defined_DM == "T2DM", "T2DM", "Non_T2DM"), levels = c("Non_T2DM", "T2DM"))
  )

# ── Mesothelial pseudobulk -> IL6, CXCL8 ───────────────────────
meso_raw <- read_tsv(meso_pb_path, show_col_types = FALSE)
gene_col <- if ("Gene" %in% colnames(meso_raw)) "Gene" else colnames(meso_raw)[1]

meso_mat <- meso_raw %>%
  mutate(across(-all_of(gene_col), as.numeric)) %>%
  column_to_rownames(gene_col) %>%
  as.matrix()

TARGET_GENES <- intersect(TARGET_GENES, rownames(meso_mat))
cat("[INFO] Target genes found:", paste(TARGET_GENES, collapse = ", "), "\n")

meso_expr_df <- as.data.frame(t(meso_mat[TARGET_GENES, , drop = FALSE])) %>%
  rownames_to_column(DONOR_COL)

# ── Mesothelial IR CASCADE module score ────────────────────────
meso_ms_df <- read_tsv(meso_ms_path, show_col_types = FALSE) %>%
  dplyr::select(all_of(c(MS_DONOR_COL, MESO_MS_COL))) %>%
  rename(!!DONOR_COL := all_of(MS_DONOR_COL), Meso_IR = all_of(MESO_MS_COL)) %>%
  distinct(.data[[DONOR_COL]], .keep_all = TRUE)

# ── align donors ───────────────────────────────────────────────
donors_common <- Reduce(intersect, list(
  meso_expr_df[[DONOR_COL]],
  meso_ms_df[[DONOR_COL]],
  cov_df[[DONOR_COL]]
)) %>% sort()

cat("[INFO] Common donors after all filters:", length(donors_common), "\n")

expr_sub <- meso_expr_df %>%
  filter(.data[[DONOR_COL]] %in% donors_common) %>%
  arrange(match(.data[[DONOR_COL]], donors_common))
ms_sub <- meso_ms_df %>%
  filter(.data[[DONOR_COL]] %in% donors_common) %>%
  arrange(match(.data[[DONOR_COL]], donors_common))
cov_sub <- cov_df %>%
  filter(.data[[DONOR_COL]] %in% donors_common) %>%
  arrange(match(.data[[DONOR_COL]], donors_common)) %>%
  dplyr::select(all_of(COV_COLS))

# ── partial correlation ────────────────────────────────────────
pcor_residual <- function(x, y, cov_df) {
  df <- data.frame(x = as.numeric(x), y = as.numeric(y), cov_df)
  df <- df[complete.cases(df), , drop = FALSE]
  n  <- nrow(df)
  if (n < MIN_N) {
    return(list(
      r = NA_real_, p = NA_real_, n = n, beta_raw = NA_real_,
      df_resid = NA_integer_, x_raw = numeric(0), y_raw = numeric(0),
      x_resid = numeric(0), y_resid = numeric(0)
    ))
  }

  X <- stats::model.matrix(~ ., data = df[, -(1:2), drop = FALSE])
  e_x <- resid(lm(df$x ~ X - 1))
  e_y <- resid(lm(df$y ~ X - 1))
  df_resid <- n - ncol(X)
  if (df_resid < 2) {
    return(list(
      r = NA_real_, p = NA_real_, n = n, beta_raw = NA_real_,
      df_resid = df_resid, x_raw = df$x, y_raw = df$y,
      x_resid = e_x, y_resid = e_y
    ))
  }

  r <- unname(cor(e_x, e_y))
  tstat <- r * sqrt(df_resid / max(1e-12, 1 - r^2))
  p <- 2 * stats::pt(-abs(tstat), df = df_resid)
  beta_raw <- unname(coef(lm(df$y ~ df$x))["df$x"])

  list(
    r = r, p = p, n = n, beta_raw = beta_raw, df_resid = df_resid,
    x_raw = df$x, y_raw = df$y, x_resid = e_x, y_resid = e_y
  )
}

res_list <- lapply(TARGET_GENES, function(g) {
  keep <- !is.na(expr_sub[[g]]) & expr_sub[[g]] != 0
  tmp <- pcor_residual(
    x = expr_sub[[g]][keep],
    y = ms_sub$Meso_IR[keep],
    cov_df = cov_sub[keep, ]
  )

  list(
    stat = data.frame(
      meso_gene = g, ir_col = MESO_MS_COL, r = tmp$r,
      beta_raw = tmp$beta_raw, p = tmp$p, n = tmp$n,
      df_resid = tmp$df_resid, stringsAsFactors = FALSE
    ),
    raw_data = data.frame(donor_id = expr_sub[[DONOR_COL]][keep], x_raw = tmp$x_raw, y_raw = tmp$y_raw),
    residual_data = data.frame(donor_id = expr_sub[[DONOR_COL]][keep], x_resid = tmp$x_resid, y_resid = tmp$y_resid)
  )
})

names(res_list) <- TARGET_GENES

res <- bind_rows(lapply(res_list, `[[`, "stat"))
raw_data <- lapply(res_list, `[[`, "raw_data")
residual_data <- lapply(res_list, `[[`, "residual_data")

cat("\n--- results ---\n")
print(res)

# ── save ───────────────────────────────────────────────────────
dir.create(dirname(OUT_RDS), recursive = TRUE, showWarnings = FALSE)

obj <- list(
  generated_at = as.character(Sys.time()),
  target_genes = TARGET_GENES,
  ir_col = MESO_MS_COL,
  donors_meso_pass = donors_meso_pass,
  donors_dp_pass = donors_dp_pass,
  donors_cell_pass = donors_cell_pass,
  donors_common = donors_common,
  covariates = COV_COLS,
  min_meso_cells = MIN_CELLS,
  min_dp_cells = MIN_CELLS,
  min_n = MIN_N,
  expr_sub = expr_sub,
  ms_sub = ms_sub,
  cov_sub = cov_sub,
  res = res,
  raw_data = raw_data,
  residual_data = residual_data
)

saveRDS(obj, OUT_RDS)
cat(sprintf("\n[SAVED] %s\n", OUT_RDS))
