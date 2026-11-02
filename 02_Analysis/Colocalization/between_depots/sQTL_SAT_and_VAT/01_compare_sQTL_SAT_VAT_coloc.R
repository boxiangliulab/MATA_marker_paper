#!/usr/bin/env Rscript
# sQTL colocalization SAT vs VAT.
#
# Decide, per cell type, whether the sQTL signals detected in both depots
# are driven by the same causal variant. Splice sites that are significant
# (global FDR < 0.05) in both SAT and VAT are tested with coloc.abf; sites
# with PP.H4 > 0.8 are called shared, the remaining significant sites are
# reported as depot-specific.
#
# coloc is run in the "eQTL mode" that uses p-values and MAF only (no sdY).
# This is more robust here than the beta/varbeta parameterisation, which
# breaks down for splice sites with standard errors of zero.

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(readr)
  library(purrr)
  library(tidyr)
  library(coloc)
  library(data.table)
  library(ggplot2)
  library(gridExtra)
})

# ====== Configuration ======
# Paths are relative to this script's directory.
results_dir <- "../../data/05_sQTL_mapping/results"
pheno_dir   <- "../../data/04_phenotype"
# Cohort MAF per population
maf_file    <- "../../data/reference/MAF_3pop_merged.tsv"
out_dir     <- "../../data/07_colocalization/sQTL_SAT_and_VAT"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Covariates used in the mapping stage: 5 splicing PCs, 10 (SAT) / 6 (VAT)
# genotype PCs; these appear in the result file names.
PC_SPLICING <- 5
G_SAT <- 10
G_VAT <- 6

# Number of donors (sample size passed to coloc)
N_DONORS <- 296

# Colocalization thresholds
PP4_THRESHOLD <- 0.8   # posterior probability of a shared causal variant
MIN_VARIANTS  <- 50    # minimum number of variants shared by the two depots

# ====== Load MAF data ======
cat("Loading MAF data...\n")
MAF <- fread(maf_file)
MAF$variant_id <- gsub("^chr", "", MAF$variant_id)
MAF[, MAF_Chinese := Chinese]
MAF[, MAF_Indian  := Indian]
MAF[, MAF_Malay   := Malay]
MAF[, MAF_mean := rowMeans(.SD, na.rm = TRUE),
    .SDcols = c("MAF_Chinese", "MAF_Indian", "MAF_Malay")]
MAF <- MAF[, .(variant_id, rsID, variant, MAF_Chinese, MAF_Indian, MAF_Malay,
               MAF_mean)]

# ====== Helper functions ======

# Number of phenotype PCs used for a depot
g_for_depot <- function(depot) if (depot == "SAT") G_SAT else G_VAT

# Read the splice sites that are significant in one depot (global FDR < 0.05)
read_sqtl_fdr05 <- function(file_path, depot) {
  if (!file.exists(file_path)) return(NULL)

  df <- read_tsv(file_path, show_col_types = FALSE, progress = FALSE) %>%
    mutate(depot = depot)

  # SNP IDs are chr:pos:ref:alt
  snp_parts <- str_split_fixed(df$SNP, ":", 4)
  df$chr <- snp_parts[, 1]
  df$pos <- as.integer(snp_parts[, 2])
  df$ref <- snp_parts[, 3]
  df$alt <- snp_parts[, 4]

  df %>%
    filter(!is.na(chr), !is.na(pos), !is.na(ref), !is.na(alt),
           chr != "", ref != "", alt != "")
}

# Read all cis variants tested for one splice site.
# The full result files are large, so the site is pulled out with grep instead
# of loading the whole file into memory.
read_all_variants_for_site <- function(cell_type, depot, splice_site, results_dir) {
  result_file <- file.path(
    results_dir,
    paste0(cell_type, ".", depot, ".p", PC_SPLICING, "_g", g_for_depot(depot),
           ".result"))

  if (!file.exists(result_file)) {
    warning("Result file not found: ", result_file)
    return(NULL)
  }

  # Note: "+" is not special in basic grep regular expressions, so splice-site
  # IDs need no escaping here.
  grep_cmd <- paste0("grep '^", splice_site, "\t' ", result_file)

  df <- tryCatch({
    read_tsv(pipe(grep_cmd),
             col_names = c("splice_site", "SNP", "pvalue", "beta", "se"),
             col_types = cols(splice_site = col_character(),
                              SNP = col_character(),
                              pvalue = col_double(),
                              beta = col_double(),
                              se = col_double()),
             show_col_types = FALSE, progress = FALSE)
  }, error = function(e) NULL)

  if (is.null(df) || nrow(df) == 0) return(NULL)

  snp_parts <- str_split_fixed(df$SNP, ":", 4)
  df$chr <- snp_parts[, 1]
  df$pos <- as.integer(snp_parts[, 2])
  df$ref <- snp_parts[, 3]
  df$alt <- snp_parts[, 4]

  df <- df %>%
    filter(!is.na(chr), !is.na(pos), !is.na(ref), !is.na(alt),
           chr != "", ref != "", alt != "")

  if (nrow(df) == 0) return(NULL)
  df
}

# Standard deviation of the PSI values of a splice site. Kept for reference;
# the coloc call below uses the p-value parameterisation and does not need sdY.
pheno_cache <- new.env()

get_sdY <- function(cell_type, depot, splice_site) {
  cache_key <- paste0(cell_type, ".", depot)

  if (!exists(cache_key, envir = pheno_cache)) {
    pheno_file <- file.path(pheno_dir, "phenotype_output", cell_type, depot,
                            "phenotype_output.noFilter_rename.prop")
    if (!file.exists(pheno_file)) {
      warning("Phenotype file not found: ", pheno_file)
      return(NA)
    }
    cat("    Loading phenotype file for", depot, "...\n")
    # Skip the sample-name line; the first column holds the splice-site ID.
    pheno_data <- fread(pheno_file, header = FALSE, skip = 1)
    setnames(pheno_data, 1, "splice_site")
    assign(cache_key, pheno_data, envir = pheno_cache)
  }

  pheno_data <- get(cache_key, envir = pheno_cache)
  site_row <- pheno_data[pheno_data$splice_site == splice_site, ]

  if (nrow(site_row) == 0) {
    warning("Splice site not found in phenotype file: ", splice_site)
    return(NA)
  }

  psi_values <- as.numeric(site_row[1, -1, with = FALSE])
  sd(psi_values, na.rm = TRUE)
}

# Colocalize the SAT and VAT signals of one splice site
run_coloc_for_site <- function(site_id, cell_type, results_dir, progress_str = "") {

  cat(progress_str, "Coloc:", site_id, "\n")

  sat_full <- read_all_variants_for_site(cell_type, "SAT", site_id, results_dir)
  vat_full <- read_all_variants_for_site(cell_type, "VAT", site_id, results_dir)

  if (is.null(sat_full) || is.null(vat_full)) {
    return(list(PP.H4.abf = NA, nsnps = 0, status = "missing_data"))
  }
  if (nrow(sat_full) == 0 || nrow(vat_full) == 0) {
    return(list(PP.H4.abf = NA, nsnps = 0, status = "no_variants"))
  }

  # Variants tested in both depots
  merged <- inner_join(
    sat_full %>% select(chr, pos, ref, alt, SNP, pvalue, beta, se),
    vat_full %>% select(chr, pos, ref, alt, SNP, pvalue, beta, se),
    by = c("chr", "pos", "ref", "alt"),
    suffix = c("_SAT", "_VAT")
  )

  if (nrow(merged) < MIN_VARIANTS) {
    return(list(PP.H4.abf = NA, nsnps = nrow(merged),
                status = "insufficient_variants"))
  }

  # Attach the cohort MAF (matched on chr:pos)
  merged$variant_id <- paste0(gsub("chr", "", merged$chr), ":", merged$pos)
  merged_dt <- as.data.table(merged)
  merged_dt <- merge(merged_dt, MAF, by = "variant_id", all.x = TRUE)
  merged_dt <- merged_dt[!is.na(MAF_Chinese) & MAF_Chinese > 0 & MAF_Chinese < 1]

  if (nrow(merged_dt) < MIN_VARIANTS) {
    return(list(PP.H4.abf = NA, nsnps = nrow(merged_dt),
                status = "insufficient_variants_after_maf"))
  }

  # Drop incomplete records
  merged_dt <- merged_dt[!is.na(pvalue_SAT) & !is.na(pvalue_VAT) &
                         !is.na(beta_SAT)  & !is.na(beta_VAT) &
                         !is.na(se_SAT)    & !is.na(se_VAT) &
                         se_SAT > 0 & se_VAT > 0]

  # One record per variant; keep the most significant when duplicated
  merged_dt <- merged_dt[order(pvalue_SAT)][!duplicated(SNP_SAT)]
  merged_dt <- merged_dt[order(pvalue_VAT)][!duplicated(SNP_VAT)]

  if (nrow(merged_dt) < MIN_VARIANTS) {
    return(list(PP.H4.abf = NA, nsnps = nrow(merged_dt),
                status = "insufficient_variants_after_filter"))
  }

  # coloc datasets: p-values + MAF (no sdY), quantitative traits
  dataset1 <- list(pvalues = merged_dt$pvalue_SAT, N = N_DONORS,
                   type = "quant", snp = merged_dt$SNP_SAT)
  dataset2 <- list(pvalues = merged_dt$pvalue_VAT, N = N_DONORS,
                   type = "quant", snp = merged_dt$SNP_VAT)

  tryCatch({
    result <- coloc.abf(dataset1, dataset2, MAF = merged_dt$MAF_Chinese)
    list(PP.H0.abf = result$summary["PP.H0.abf"],
         PP.H1.abf = result$summary["PP.H1.abf"],
         PP.H2.abf = result$summary["PP.H2.abf"],
         PP.H3.abf = result$summary["PP.H3.abf"],
         PP.H4.abf = result$summary["PP.H4.abf"],
         nsnps = result$summary["nsnps"],
         status = "success",
         merged_data = merged_dt)      # kept for the locus plot
  }, error = function(e) {
    cat("    ERROR:", e$message, "\n")
    list(PP.H4.abf = NA, nsnps = nrow(merged_dt), status = "coloc_error")
  })
}

# Locus plot of one colocalized site: -log10(p) in SAT and VAT side by side
plot_coloc_site <- function(site_id, coloc_result, plot_dir, rank, category) {
  if (is.null(coloc_result$merged_data) ||
      nrow(coloc_result$merged_data) < 10) {
    return(NULL)
  }

  merged_dt <- coloc_result$merged_data
  merged_dt[, logp_SAT := -log10(pvalue_SAT)]
  merged_dt[, logp_VAT := -log10(pvalue_VAT)]

  lead_sat <- merged_dt[which.min(pvalue_SAT)]
  lead_vat <- merged_dt[which.min(pvalue_VAT)]

  # Agreement of the effect sizes between depots
  cor_beta <- cor(merged_dt$beta_SAT, merged_dt$beta_VAT)
  merged_dt[, same_direction := sign(beta_SAT) == sign(beta_VAT)]
  pct_concordant <- 100 * sum(merged_dt$same_direction) / nrow(merged_dt)

  plot_data <- merged_dt[, .(pos, logp_SAT, logp_VAT)]
  is_lead_SAT <- merged_dt$pos == lead_sat$pos
  is_lead_VAT <- merged_dt$pos == lead_vat$pos
  plot_data[, point_type := "Other"]
  plot_data[is_lead_SAT, point_type := "Lead_SAT"]
  plot_data[is_lead_VAT, point_type := "Lead_VAT"]
  plot_data[is_lead_SAT & is_lead_VAT, point_type := "Lead_Both"]

  chrom <- unique(merged_dt$chr)[1]

  point_colors <- c("Other" = "gray50", "Lead_SAT" = "red",
                    "Lead_VAT" = "blue", "Lead_Both" = "purple")
  point_sizes <- c("Other" = 1, "Lead_SAT" = 3,
                   "Lead_VAT" = 2, "Lead_Both" = 4)

  p1 <- ggplot(plot_data, aes(x = pos / 1e6, y = logp_SAT)) +
    geom_point(aes(color = point_type, size = point_type), alpha = 0.7) +
    scale_color_manual(values = point_colors) +
    scale_size_manual(values = point_sizes) +
    labs(title = paste0(site_id, " - SAT"),
         subtitle = paste0(chrom, " | Effect correlation: ",
                           round(cor_beta, 3)),
         x = "Position (Mb)", y = "-log10(P-value)") +
    theme_bw(base_size = 11) +
    theme(legend.position = "none", panel.grid.minor = element_blank())

  p2 <- ggplot(plot_data, aes(x = pos / 1e6, y = logp_VAT)) +
    geom_point(aes(color = point_type, size = point_type), alpha = 0.7) +
    scale_color_manual(values = point_colors, name = "Variant type") +
    scale_size_manual(values = point_sizes, name = "Variant type") +
    labs(title = paste0(site_id, " - VAT"),
         subtitle = paste0("Sign concordance: ",
                           round(pct_concordant, 1), "%"),
         x = "Position (Mb)", y = "-log10(P-value)") +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())

  site_name <- gsub(":", "_", site_id)
  output_file <- file.path(plot_dir,
                           paste0("coloc_example_", rank, "_", site_name, ".pdf"))

  pdf(output_file, width = 10, height = 8)
  grid.arrange(p1, p2, ncol = 1,
               top = paste0("Example ", rank, " - ", category, ": ",
                            "PP.H0 = ", round(coloc_result$PP.H0.abf, 4), " | ",
                            "PP.H3 = ", round(coloc_result$PP.H3.abf, 4), " | ",
                            "PP.H4 = ", round(coloc_result$PP.H4.abf, 4)))
  dev.off()

  TRUE
}

# ====== Analysis of one cell type ======

analyse_cell_type <- function(CELL_TYPE) {


  cat("\n========================================\n")
  cat("Cell type:", CELL_TYPE, "\n")
  cat("========================================\n")

  sat_file <- file.path(results_dir,
                        paste0(CELL_TYPE, ".SAT.p", PC_SPLICING, "_g", G_SAT,
                               ".result.sQTL.fdr05"))
  vat_file <- file.path(results_dir,
                        paste0(CELL_TYPE, ".VAT.p", PC_SPLICING, "_g", G_VAT,
                               ".result.sQTL.fdr05"))

  if (!file.exists(sat_file)) stop("SAT file not found: ", sat_file)
  if (!file.exists(vat_file)) stop("VAT file not found: ", vat_file)

  sat_sqtl <- read_sqtl_fdr05(sat_file, "SAT")
  vat_sqtl <- read_sqtl_fdr05(vat_file, "VAT")

  if (is.null(sat_sqtl) || is.null(vat_sqtl)) stop("Failed to read sQTL data")

  sat_sites <- unique(sat_sqtl$splice_site)
  vat_sites <- unique(vat_sqtl$splice_site)

  cat("  SAT sQTL sites (FDR < 0.05):", length(sat_sites), "\n")
  cat("  VAT sQTL sites (FDR < 0.05):", length(vat_sites), "\n")

  # Sites significant in both depots are the colocalization candidates
  shared_sites <- intersect(sat_sites, vat_sites)
  cat("  Shared sites (candidates for colocalization):", length(shared_sites), "\n")

  ct_out <- file.path(out_dir, CELL_TYPE)
  dir.create(ct_out, recursive = TRUE, showWarnings = FALSE)

  if (length(shared_sites) == 0) {
    cat("  No shared sites, saving empty results\n")

    write_tsv(tibble(splice_site = character(0)),
              file.path(ct_out, "colocalized.tsv"))
    write_tsv(tibble(splice_site = sat_sites), file.path(ct_out, "SAT_specific.tsv"))
    write_tsv(tibble(splice_site = vat_sites), file.path(ct_out, "VAT_specific.tsv"))

    summary_df <- tibble(cell_type = CELL_TYPE,
                         n_SAT = length(sat_sites),
                         n_VAT = length(vat_sites),
                         n_shared_candidates = 0,
                         n_colocalized = 0,
                         n_SAT_specific = length(sat_sites),
                         n_VAT_specific = length(vat_sites))
    write_tsv(summary_df, file.path(ct_out, "summary.tsv"))

    return(summary_df)
  }

  plot_dir <- file.path(ct_out, "plots")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

  cat("  Running colocalization...\n")
  n_shared <- length(shared_sites)

  coloc_results <- map_dfr(seq_along(shared_sites), function(i) {
    site <- shared_sites[i]
    progress <- sprintf("    [%d/%d] ", i, n_shared)
    res <- run_coloc_for_site(site, CELL_TYPE, results_dir, progress)

    if (res$status == "success" && !is.null(res$merged_data)) {
      tryCatch({
        plot_coloc_site(site, res, plot_dir, i, "coloc")
      }, error = function(e) {
        cat("    Warning: Failed to plot:", e$message, "\n")
      })
    }

    tibble(splice_site = site,
           PP.H0.abf = res$PP.H0.abf %||% NA,
           PP.H1.abf = res$PP.H1.abf %||% NA,
           PP.H2.abf = res$PP.H2.abf %||% NA,
           PP.H3.abf = res$PP.H3.abf %||% NA,
           PP.H4.abf = res$PP.H4.abf %||% NA,
           nsnps = res$nsnps,
           status = res$status)
  })

  # ====== Classify the splice sites ======
  colocalized_sites <- coloc_results %>%
    filter(!is.na(PP.H4.abf), PP.H4.abf > PP4_THRESHOLD) %>%
    pull(splice_site)

  not_colocalized <- setdiff(shared_sites, colocalized_sites)
  sat_specific <- setdiff(sat_sites, c(colocalized_sites, not_colocalized))
  vat_specific <- setdiff(vat_sites, c(colocalized_sites, not_colocalized))

  cat("  Colocalized (PP.H4 >", PP4_THRESHOLD, "):", length(colocalized_sites), "\n")
  cat("  SAT-specific:", length(sat_specific), "\n")
  cat("  VAT-specific:", length(vat_specific), "\n")

  write_tsv(tibble(splice_site = colocalized_sites),
            file.path(ct_out, "colocalized.tsv"))
  write_tsv(tibble(splice_site = sat_specific), file.path(ct_out, "SAT_specific.tsv"))
  write_tsv(tibble(splice_site = vat_specific), file.path(ct_out, "VAT_specific.tsv"))
  write_tsv(coloc_results, file.path(ct_out, "coloc_results_full.tsv"))

  summary_df <- tibble(cell_type = CELL_TYPE,
                       n_SAT = length(sat_sites),
                       n_VAT = length(vat_sites),
                       n_shared_candidates = length(shared_sites),
                       n_colocalized = length(colocalized_sites),
                       n_SAT_specific = length(sat_specific),
                       n_VAT_specific = length(vat_specific))
  write_tsv(summary_df, file.path(ct_out, "summary.tsv"))

  cat("\nResults saved to:", ct_out, "\n")
  cat("Plots saved to:", plot_dir, "\n")
  cat("========================================\n")

  return(summary_df)
}

# ====== Entry point ======
# With one or more cell types on the command line only those are processed;
# with no argument every cell type tested in both depots is processed in turn
# and the per-cell-type summaries are combined.

ref_dir <- "../../reference"

celltypes_both_depots <- function() {
  read_list <- function(depot) {
    f <- file.path(ref_dir, paste0("celltypes.QTL.", depot, ".txt"))
    if (!file.exists(f)) stop("Cell-type list not found: ", f)
    read.table(f, sep = "\t", stringsAsFactors = FALSE)[[1]]
  }
  intersect(read_list("SAT"), read_list("VAT"))
}

args <- commandArgs(trailingOnly = TRUE)
run_all <- length(args) == 0
cell_types <- if (run_all) celltypes_both_depots() else args

summaries <- lapply(cell_types, analyse_cell_type)

if (run_all) {
  summary_all <- bind_rows(summaries)
  summary_file <- file.path(out_dir, "summary_all_celltypes.tsv")
  write_tsv(summary_all, summary_file)
  cat("\nCombined summary written to:", summary_file, "\n")
}
