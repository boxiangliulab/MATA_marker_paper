#!/usr/bin/env Rscript
# Purpose: Test mesothelial fibrosis-gene associations using continuous DP-ATM
# fraction rather than binary DP-ATM status.
# Inputs: SECTION4_MESOTHELIAL_PSEUDOBULK, SECTION4_VAT_ATM_COUNTS, and
# SECTION4_METADATA.
# Output: Model estimates in TSV format.
# Method: Linear models for a 10-percentage-point increase in DP fraction,
# adjusted for donor covariates, with and without total VAT cell count.
# Usage: Rscript mesothelial_fibrosis_continuous_dp_fraction.R
suppressPackageStartupMessages({ library(data.table); library(dplyr)
                                 library(tidyr); library(broom) })
out_dir <- file.path(Sys.getenv("SECTION4_OUTPUT"), "sensitivity")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
GENES <- c("COL1A1", "CCN2", "SERPINE1", "TGFBI", "FN1", "TNC")

expr <- fread(Sys.getenv("SECTION4_MESOTHELIAL_PSEUDOBULK"))
gcol <- names(expr)[1]
e <- melt(expr[get(gcol) %in% GENES], id.vars = gcol,
          variable.name = "donor_id", value.name = "expr")
setnames(e, gcol, "gene"); e[, donor_id := as.character(donor_id)]
e <- dcast(e, donor_id ~ gene, value.var = "expr")

cnt <- fread(Sys.getenv("SECTION4_VAT_ATM_COUNTS"))
cov <- unique(fread(Sys.getenv("SECTION4_METADATA"),
                    select = c("donor_id","Age","Gender","Race","BMI","Self_defined_DM")),
              by = "donor_id")
num <- function(v) suppressWarnings(as.numeric(ifelse(as.character(v)=="x", NA, as.character(v))))
cov[, `:=`(Age = num(Age), BMI = num(BMI), Gender = factor(Gender), Race = factor(Race),
           # Match the binary disease coding used in the primary analysis.
           DM = factor(ifelse(Self_defined_DM == "T2DM", "T2DM", "Non_T2DM")))]

d <- Reduce(function(a, b) merge(a, b, by = "donor_id"), list(e, cnt, cov))
d <- d[n_ATM >= 20]                      # a fraction needs a usable denominator
d[, `:=`(DP_frac10 = 10 * n_DP / n_ATM,  # per 10 percentage points
         log10_cells = log10(n_VAT_cells))]
d <- d[complete.cases(d[, .(Age, Gender, Race, BMI, DM)])]
cat("donors used:", nrow(d), "\n")
cat("DP fraction of ATMs: median", round(median(d$n_DP/d$n_ATM), 3),
    " range", paste(round(range(d$n_DP/d$n_ATM), 3), collapse = " to "), "\n")
cat("cor(DP fraction, VAT cells captured) Spearman r =",
    round(cor(d$DP_frac10, d$n_VAT_cells, method = "spearman"), 3), "\n\n")

RHS <- c(cont      = "DP_frac10 + Age + Gender + Race + BMI + DM",
         cont_cells= "DP_frac10 + Age + Gender + Race + BMI + DM + n_VAT_cells")
run <- function(g, rhs, nm) {
  dat <- d %>% select(all_of(c(g, "DP_frac10", "Age", "Gender", "Race", "BMI",
                               "DM", "n_VAT_cells"))) %>% drop_na()
  tidy(lm(as.formula(paste(g, "~", rhs)), data = dat), conf.int = TRUE) %>%
    filter(term == "DP_frac10") %>%
    transmute(gene = g, model = nm, n = nrow(dat), beta = estimate,
              se = std.error, p = p.value)
}
res <- rbindlist(lapply(names(RHS), function(nm)
  rbindlist(lapply(GENES, run, rhs = RHS[[nm]], nm = nm))))
res[, padj := p.adjust(p, "BH"), by = model]
res[, sig := cut(padj, c(-Inf,.001,.01,.05,Inf), c("***","**","*","ns"))]
res[, gene := factor(gene, levels = GENES)]
for (nm in names(RHS)) {
  cat("---", nm, ":", RHS[[nm]], "\n")
  print(res[model == nm][order(gene),
            .(gene, n, beta_per10pct = round(beta,3), p = signif(p,3),
              padj = signif(padj,3), sig)]); cat("\n")
}
fwrite(res, file.path(out_dir, "mesothelial_fibrosis_continuous_dp_fraction.tsv"), sep = "\t")
