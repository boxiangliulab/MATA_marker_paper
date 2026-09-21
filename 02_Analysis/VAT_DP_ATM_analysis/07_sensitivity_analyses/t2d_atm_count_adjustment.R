#!/usr/bin/env Rscript
# Purpose: Test whether the association between DP-ATM status and T2D is robust
# to VAT ATM abundance and donor covariates.
# Inputs: SECTION4_T2D_STATUS_RDS and SECTION4_VAT_ATM_COUNTS.
# Outputs: Logistic-regression estimates in TSV and a complete analysis RDS.
# Method: Logistic regression with T2D as the outcome and alternative ATM count,
# VAT cell count, and donor-covariate adjustments.
# Usage: Rscript t2d_atm_count_adjustment.R
suppressPackageStartupMessages({ library(data.table); library(dplyr); library(broom) })
out_dir <- file.path(Sys.getenv("SECTION4_OUTPUT"), "sensitivity")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

o  <- readRDS(Sys.getenv("SECTION4_T2D_STATUS_RDS"))
d  <- as.data.table(o$donor_df)[!is.na(group_DP) & !is.na(DM_binary)]
cnt <- fread(Sys.getenv("SECTION4_VAT_ATM_COUNTS"))
d <- merge(d, cnt[, .(donor_id, n_ATM, n_VAT_cells, n_DP_freeze = n_DP)], by = "donor_id")
cat("donors:", nrow(d), " (primary analysis n =", o$chisq_res_dp$n_total, ")\n")
cat("n_DP matches the freeze recount for", sum(d$n_DP == d$n_DP_freeze), "of", nrow(d), "\n\n")

num <- function(v) suppressWarnings(as.numeric(ifelse(as.character(v) == "x", NA, as.character(v))))
d[, `:=`(Age = num(Age), BMI = num(BMI), Gender = factor(Gender), Race = factor(Race),
         group_DP = factor(group_DP, levels = c("DP absent", "DP present")),
         T2D = as.integer(DM_binary == "T2DM"),
         ATM100 = n_ATM / 100, cells1k = n_VAT_cells / 1000)]

cat("primary chi-square: p =", signif(o$chisq_res_dp$p_value, 3), "\n")
print(o$contingency_dp)
cat("\nATM count by DP group: median",
    paste(d[, as.numeric(median(n_ATM)), by = group_DP][order(group_DP)]$V1, collapse = " vs "),
    " Wilcoxon p =", signif(wilcox.test(n_ATM ~ group_DP, data = d)$p.value, 3), "\n")
cat("ATM fraction        : median",
    paste(round(d[, median(n_ATM / n_VAT_cells), by = group_DP][order(group_DP)]$V1, 3), collapse = " vs "),
    " Wilcoxon p =", signif(wilcox.test(I(n_ATM/n_VAT_cells) ~ group_DP, data = d)$p.value, 3), "\n\n")

MOD <- list(
  L1_crude       = "group_DP",
  L2_ATM         = "group_DP + ATM100",
  L3_cells       = "group_DP + cells1k",
  L4_cov         = "group_DP + Age + Gender + Race + BMI",
  L5_cov_ATM     = "group_DP + Age + Gender + Race + BMI + ATM100",
  L6_cov_cells   = "group_DP + Age + Gender + Race + BMI + cells1k")

res <- rbindlist(lapply(names(MOD), function(nm) {
  dat <- d %>% select(T2D, group_DP, Age, Gender, Race, BMI, ATM100, cells1k) %>% tidyr::drop_na()
  fit <- glm(as.formula(paste("T2D ~", MOD[[nm]])), data = dat, family = binomial)
  tidy(fit, conf.int = TRUE) %>% filter(term == "group_DPDP present") %>%
    transmute(model = nm, formula = MOD[[nm]], n = nrow(dat),
              OR = exp(estimate), lo = exp(conf.low), hi = exp(conf.high), p = p.value)
}))
res[, sig := cut(p, c(-Inf, .001, .01, .05, Inf), c("***", "**", "*", "ns"))]
print(res[, .(model, n, OR = round(OR, 2), CI = sprintf("%.2f-%.2f", lo, hi),
              p = signif(p, 3), sig)])
saveRDS(list(generated_at = as.character(Sys.time()), donors = d,
             chisq_primary = o$chisq_res_dp, result = res),
        file.path(out_dir, "t2d_atm_count_adjustment.rds"))
fwrite(res, file.path(out_dir, "t2d_atm_count_adjustment.tsv"), sep = "\t")
cat("\n[saved] t2d_atm_count_adjustment.{rds,tsv}\n")
