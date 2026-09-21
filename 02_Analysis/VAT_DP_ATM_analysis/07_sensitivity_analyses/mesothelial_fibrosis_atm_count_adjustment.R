#!/usr/bin/env Rscript
# Purpose: Test whether associations between DP-ATM status and mesothelial
# fibrosis-gene expression are robust to VAT ATM and cell-recovery measures.
# Inputs: SECTION4_MESOTHELIAL_FIBROSIS_RDS and SECTION4_VAT_ATM_COUNTS.
# Outputs: model estimates in TSV and a complete analysis object in RDS format.
# Method: Linear models adjusted for donor covariates, with alternative ATM
# count, VAT cell count, ATM fraction, and log10 ATM count adjustments.
# Usage: Rscript mesothelial_fibrosis_atm_count_adjustment.R
suppressPackageStartupMessages({ library(data.table); library(dplyr)
                                 library(tidyr); library(broom) })
out_dir <- file.path(Sys.getenv("SECTION4_OUTPUT"), "sensitivity")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
o <- readRDS(Sys.getenv("SECTION4_MESOTHELIAL_FIBROSIS_RDS"))
GENES <- o$genes
cnt <- fread(Sys.getenv("SECTION4_VAT_ATM_COUNTS"))

d <- as.data.table(o$df_dp)[, .(donor_id, COL1A1, CCN2, SERPINE1, TGFBI, FN1, TNC,
                                Age, Gender, Race, BMI, Self_defined_DM,
                                n_DP_orig = n_DP, group_DP)]
d <- merge(d, cnt, by = "donor_id", all.x = TRUE)
cat("donors:", nrow(d), " matched to VAT counts:", sum(!is.na(d$n_ATM)), "\n")
cat("n_DP in the input RDS vs recomputed counts: identical for",
    sum(d$n_DP_orig == d$n_DP, na.rm = TRUE), "of", sum(!is.na(d$n_DP)), "donors\n\n")

d[, `:=`(Gender = factor(Gender), Race = factor(Race),
         Self_defined_DM = factor(Self_defined_DM),
         log10_ATM = log10(n_ATM), frac_ATM = n_ATM / n_VAT_cells,
         log10_cells = log10(n_VAT_cells))]

# Examine the relation between cell-recovery measures and DP status.
dd <- d[!is.na(Self_defined_DM) & Self_defined_DM != "Uncertain"]
dd[, Self_defined_DM := droplevels(Self_defined_DM)]
cat("ATM count by DP group\n"); print(dd[, .(n = .N, median_ATM = as.numeric(median(n_ATM)),
      median_frac = round(median(frac_ATM), 3)), by = group_DP])
cat("  Wilcoxon n_ATM  DP_present vs DP_absent  p =",
    signif(wilcox.test(n_ATM ~ group_DP, data = dd)$p.value, 3), "\n\n")

BASE <- "group_DP + Age + Gender + Race + BMI + Self_defined_DM"
RHS <- c(
  M2        = BASE,
  M3_ATM    = paste(BASE, "+ n_ATM"),
  M3_cells  = paste(BASE, "+ n_VAT_cells"),
  M3_frac   = paste(BASE, "+ frac_ATM"),
  M3_log    = paste(BASE, "+ log10_ATM"))

# Report coefficients in both log2 CPM and gene-specific standard-deviation units.
fit1 <- function(g, rhs, nm) {
  dat <- dd %>% select(all_of(c(g, "group_DP", "Age", "Gender", "Race", "BMI",
                                "Self_defined_DM", "n_ATM", "n_VAT_cells",
                                "log10_ATM", "frac_ATM"))) %>% drop_na()
  sdv <- sd(dat[[g]])
  dat[[paste0(g, "_z")]] <- dat[[g]] / sdv
  grab <- function(y) tidy(lm(as.formula(paste(y, "~", rhs)), data = dat),
                           conf.int = TRUE) %>% filter(term == "group_DPDP_present")
  a <- grab(g); b <- grab(paste0(g, "_z"))
  transmute(a, gene = g, model = nm, n = nrow(dat), beta = estimate,
            se = std.error, lo = conf.low, hi = conf.high, p = p.value,
            beta_sd = b$estimate, lo_sd = b$conf.low, hi_sd = b$conf.high)
}
res <- rbindlist(lapply(names(RHS), function(nm)
  rbindlist(lapply(GENES, fit1, rhs = RHS[[nm]], nm = nm))))
res[, gene := factor(gene, levels = GENES)]
res[, padj := p.adjust(p, "BH"), by = model]
res[, sig := cut(padj, c(-Inf, .001, .01, .05, Inf), c("***", "**", "*", "ns"))]

for (nm in names(RHS)) {
  cat("---", nm, ":", RHS[[nm]], "\n")
  print(res[model == nm, .(gene, n, beta = round(beta, 3), p = signif(p, 3),
                           padj = signif(padj, 3), sig)])
  cat("\n")
}
saveRDS(list(generated_at = as.character(Sys.time()), formulas = RHS,
             donors = dd, result = res),
        file.path(out_dir, "mesothelial_fibrosis_atm_count_adjustment.rds"))
fwrite(res, file.path(out_dir, "mesothelial_fibrosis_atm_count_adjustment.tsv"), sep = "\t")
cat("[saved] mesothelial_fibrosis_atm_count_adjustment.{rds,tsv}\n")
