#!/usr/bin/env Rscript
# Colocalization of the sQTL with GWAS signals.
#
# For one cell type and depot, test whether the cis-sQTL signals share a
# causal variant with published GWAS signals.
#
# For every GWAS trait the script
#   1. intersects the QTL loci with the genome-wide significant GWAS variants
#      (p < 1e-7 subset), drops the HLA region (chr6:25-35 Mb), and clumps the
#      remaining variants so that no two retained signals are within 500 kb;
#   2. runs coloc.abf on each retained splice site, using the full GWAS
#      summary statistics and all cis variants of that splice site, with the
#      cohort MAF; case-control and quantitative traits are handled according
#      to the trait table;
#   3. writes the posterior probabilities plus the statistics of the top QTL
#      variant, and stores the input summary statistics of every locus with
#      PP.H4 > 0.75 for plotting.
#
# Adapted from a GWAS-QTL colocalization script written by Yihan Tong
# (original version 2025-06-24).

suppressPackageStartupMessages({
  library(data.table)
  library(getopt)
  library(dplyr)
  library(coloc)
})

command <- matrix(c(
  "input_file",               "i", 1, "character",
  "output_inter_file",        "o", 1, "character",
  "output_coloc_file_prefix", "c", 1, "character",
  "qtl_sum_stats_dir",        "d", 1, "character",
  "N",                        "N", 1, "integer",
  "cell_type",                "t", 1, "character",
  "help",                     "h", 0, "logical"), byrow = TRUE, ncol = 4)
args <- getopt(command)

if (!is.null(args$help)) {
  cat("Usage: Rscript 01_coloc_sQTL_with_GWAS.R [options]\n")
  cat("Options:\n")
  cat("  -i, --input_file                QTL loci of one cell type/depot (required)\n")
  cat("  -o, --output_inter_file         Output file for the intersected loci (required)\n")
  cat("  -c, --output_coloc_file_prefix  Output prefix for the coloc results (required)\n")
  cat("  -d, --qtl_sum_stats_dir         Directory with the per-splice-site results (required)\n")
  cat("  -N, --N                         QTL sample size (required)\n")
  cat("  -t, --cell_type                 Restrict the batch run to one cell type\n")
  cat("  -h, --help                      Show this message\n")
  cat("\nWith no options at all, every cell type and depot listed in\n")
  cat("reference/celltypes.QTL.<depot>.txt is processed in turn.\n")
  q(status = 0)
}

# --- Reference data (paths are relative to this script's directory) ----------
maf_file         <- "../../data/reference/MAF_3pop_merged.tsv"
gwas_dir         <- "../../data/reference/GWAS_SumStats"
gwas_trait_table <- "../../data/reference/GWAS_trait.txt"

# Analysis parameters
HLA_CHR <- "6"; HLA_START <- 25e6; HLA_END <- 35e6   # HLA region, excluded
CLUMP_WINDOW <- 500000                               # clumping window (bp)
MIN_SNPS_FOR_COLOC <- 10                             # minimum overlap for coloc
PP4_SAVE_THRESHOLD <- 0.75                           # save sum-stats above this

########## Reference files ##########
# Cohort MAF (three populations); variant_id is chr:pos without the "chr"
MAF <- fread(maf_file)
MAF$variant_id <- gsub("^chr", "", MAF$variant_id)
MAF[, MAF_Chinese := Chinese]
MAF[, MAF_Indian  := Indian]
MAF[, MAF_Malay   := Malay]
MAF[, MAF_mean := rowMeans(.SD, na.rm = TRUE),
    .SDcols = c("MAF_Chinese", "MAF_Indian", "MAF_Malay")]
MAF <- MAF[, .(variant_id, rsID, variant, MAF_Chinese, MAF_Indian, MAF_Malay,
               MAF_mean)]

# GWAS sample sizes and study design
sample_size_gwas <- read.table(gwas_trait_table, sep = "\t", header = TRUE)
sample_size_gwas$Study <- gsub("_hg38", "", sample_size_gwas$Study)
sample_size_gwas$N <- gsub(",", "", sample_size_gwas$N) %>% as.numeric()

# Traits: one file per trait holding the genome-wide significant variants
gwas_sig_dir <- file.path(gwas_dir, "threshold_1e-7")
filelist <- list.files(gwas_sig_dir, pattern = "hg38_1e_7\\.txt$",
                       full.names = TRUE)
trait_names <- gsub("_hg38_1e_7\\.txt$", "", basename(filelist))

# ====== Colocalization of one cell type and depot ======
# input_file               candidate QTL loci (splice_site, SNP, pvalue)
# output_inter_file        where the intersected, clumped loci are written
# output_coloc_file_prefix prefix of the per-trait coloc result files
# qtl_sum_stats_dir        directory with the per-splice-site association files
# N_qtl                    QTL sample size

run_cell_type_depot <- function(input_file, output_inter_file,
                                output_coloc_file_prefix, qtl_sum_stats_dir,
                                N_qtl) {

  dir.create(dirname(output_inter_file), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(output_coloc_file_prefix), recursive = TRUE,
             showWarnings = FALSE)

  # --- QTL loci (loaded once, used for every trait) ---------------------------
  message("  - Loading QTL nominal loci...")
  nominal <- fread(input_file)
  message("    nominal rows: ", nrow(nominal))
  setnames(nominal, c("SNP", "pvalue", "splice_site"),
           c("variant_id", "p.value", "phenotype"))
  # Variant IDs are reduced to chr:pos (no "chr" prefix) to match the GWAS files
  nominal$variant_id <- gsub("chr", "", nominal$variant_id)
  nominal$variant_id <- sub("(:[^:]+){2}$", "", nominal$variant_id)
  message("    nominal after renaming rows: ", nrow(nominal))

  for (trait in trait_names) {
    message(">>> Processing trait: ", trait)

    message("  - Loading GWAS summary stats...")
    GWAS_pre <- fread(file.path(gwas_sig_dir, paste0(trait, "_hg38_1e_7.txt")))
    GWAS_raw <- fread(file.path(gwas_dir, paste0(trait, "_hg38.txt")))

    if (nrow(GWAS_pre) == 0) {
      message("  - No trait found in ", trait)
      next
    }

    GWAS_pre <- GWAS_pre[, .(variant_id = SNPID_hg38, beta = BETA, se = SE,
                             Pvalue = Pvalue)]
    GWAS_raw <- GWAS_raw[, .(variant_id = SNPID_hg38, beta = BETA, se = SE,
                             Pvalue = Pvalue)]
    message("  - GWAS_pre rows: ", nrow(GWAS_pre))

    # --- Variants shared by the QTL loci and the GWAS signals -------------------
    message("  - intersect variant_id count: ",
            length(intersect(nominal$variant_id, GWAS_pre$variant_id)))
    common_snp <- merge(nominal, GWAS_pre, by = "variant_id", all = FALSE)
    common_snp <- na.omit(common_snp)
    message("  - common_snp rows: ", nrow(common_snp))

    if (nrow(common_snp) == 0) {
      message("  - No intersecting SNPs found. Skipping.")
      next
    }

    # --- Remove the HLA region --------------------------------------------------
    message("  - Filtering out HLA region...")
    chr_pos <- do.call(rbind, strsplit(common_snp$variant_id, ":"))
    chr <- chr_pos[, 1]
    pos <- as.numeric(chr_pos[, 2])
    common_snp <- common_snp[!(chr == HLA_CHR & pos >= HLA_START & pos <= HLA_END), ]

    if (nrow(common_snp) == 0) {
      message("  - All SNPs filtered by HLA. Skipping.")
      next
    }

    # --- Clumping ---------------------------------------------------------------
    # Keep the most significant variants per chromosome, making sure that no two
    # retained variants lie within CLUMP_WINDOW of each other, so that nearby
    # variants of the same signal are not tested repeatedly.
    message("  - Clumping SNPs...")
    common_snp <- as.data.table(common_snp)
    common_snp[, c("chr", "pos") := tstrsplit(variant_id, ":", fixed = TRUE)]
    common_snp[, pos := as.numeric(pos)]
    common_snp[, keep := FALSE]

    for (current_chr in unique(common_snp$chr)) {
      dt_chr <- common_snp[chr == current_chr]
      dt_chr <- dt_chr[order(Pvalue)]
      dt_chr[, keep_tmp := TRUE]

      for (i in 1:nrow(dt_chr)) {
        if (dt_chr$keep_tmp[i]) {
          current_pos <- dt_chr$pos[i]
          dt_chr[pos > current_pos & pos - current_pos < CLUMP_WINDOW,
                 keep_tmp := FALSE]
        }
      }
      common_snp[chr == current_chr &
                   variant_id %in% dt_chr[keep_tmp == TRUE]$variant_id,
                 keep := TRUE]
    }

    common_snp_filtered <- common_snp[keep == TRUE]
    message("  - ", nrow(common_snp_filtered), " SNPs retained after clumping.")

    common_snp_filtered$trait <- trait
    # Shorten the trait name (drop the trailing year and journal) to keep the
    # output paths within the file-system limits.
    trait_short <- gsub("_[0-9]{4}_[A-Za-z_]+$", "", trait)
    output_inter_file_this_trait <- file.path(
      dirname(output_inter_file),
      paste0(trait_short, ".", basename(output_inter_file)))
    message("  - Output intersection file: ", output_inter_file_this_trait)

    write.table(common_snp_filtered[, !c("chr", "pos", "keep")],
                output_inter_file_this_trait,
                sep = "\t", quote = FALSE, row.names = FALSE)

    unique_phenotype_need_test <- unique(common_snp_filtered$phenotype)

    if (nrow(common_snp_filtered) > 0) {
      message("  - Starting coloc.abf for retained SNPs...")

      test_rs <- data.frame(matrix(NA, nrow = length(unique_phenotype_need_test),
                                   ncol = 17))
      colnames(test_rs) <- c("variant_id", "phenotype", "nsnps",
                             "PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf",
                             "PP.H4.abf", "Pvalue_gwas", "Pvalue_qtl",
                             "beta_gwas", "se_gwas", "beta_qtl",
                             "rsID", "MAF_Chinese", "MAF_Indian", "MAF_Malay")

      for (j in seq_along(unique_phenotype_need_test)) {
        phenotype <- unique_phenotype_need_test[j]

        # All cis variants tested for this splice site
        coloc_file <- file.path(qtl_sum_stats_dir, paste0(phenotype, ".result"))
        message("    - coloc_file: ", coloc_file,
                ", exists: ", file.exists(coloc_file))
        if (!file.exists(coloc_file)) next

        snp_1mb <- fread(coloc_file, sep = "\t", header = FALSE)
        setnames(snp_1mb, c("splice_site", "SNP", "pvalue", "beta", "se"))
        snp_1mb[, variant_id := gsub("chr", "", SNP)]
        snp_1mb[, variant_id := sapply(strsplit(variant_id, ":"),
                                       function(x) paste(x[1], x[2], sep = ":"))]
        setnames(snp_1mb, c("splice_site", "pvalue", "beta", "se"),
                 c("phenotype", "p.value_qtl", "BETA_qtl", "SE_qtl"))
        if (!"N" %in% colnames(snp_1mb)) snp_1mb$N_qtl <- N_qtl

        # Top QTL variant of this splice site
        snp_1mb <- snp_1mb[order(p.value_qtl)]
        top_snp <- snp_1mb[1, "variant_id"]
        test_rs[j, c("variant_id", "phenotype")] <- snp_1mb[1, .(variant_id,
                                                                 phenotype)]
        test_rs[j, "Pvalue_qtl"] <- snp_1mb[1, p.value_qtl]
        test_rs[j, "beta_qtl"]   <- snp_1mb[1, BETA_qtl]

        # Variants shared by the QTL, the MAF table and the GWAS
        input <- merge(snp_1mb, MAF, by = "variant_id")
        GWAS_tmp <- copy(GWAS_raw)
        setnames(GWAS_tmp, c("beta", "se", "Pvalue"),
                 c("beta_gwas", "se_gwas", "Pvalue_gwas"))
        input <- merge(input, GWAS_tmp, by = "variant_id")
        input <- input[order(Pvalue_gwas)]
        input[, varbeta := se_gwas^2]

        unique_input <- input[!duplicated(variant_id)]
        n_before_maf_filter <- nrow(unique_input)
        unique_input <- unique_input[!is.na(MAF_Chinese) & MAF_Chinese > 0 &
                                       MAF_Chinese < 1]
        n_after_maf_filter <- nrow(unique_input)
        # Incomplete GWAS records would make coloc.abf fail
        unique_input <- unique_input[!is.na(se_gwas) & !is.na(beta_gwas) &
                                       !is.na(Pvalue_gwas)]
        message("    - unique_input rows: ", nrow(unique_input),
                " (filtered: ", n_before_maf_filter - n_after_maf_filter,
                " by MAF, ", n_after_maf_filter - nrow(unique_input),
                " by missing GWAS)")

        test_rs[j, "nsnps"] <- nrow(unique_input)

        # GWAS statistics of the top QTL variant
        idx_g <- match(top_snp, GWAS_raw$variant_id)
        if (!is.na(idx_g)) {
          test_rs[j, "Pvalue_gwas"] <- GWAS_raw$Pvalue[idx_g]
          test_rs[j, "beta_gwas"]   <- GWAS_raw$beta[idx_g]
          test_rs[j, "se_gwas"]     <- GWAS_raw$se[idx_g]
        }
        # rsID and population MAF of the top QTL variant
        idx_m <- match(top_snp, unique_input$variant_id)
        if (!is.na(idx_m)) {
          test_rs[j, "rsID"]        <- unique_input$rsID[idx_m]
          test_rs[j, "MAF_Chinese"] <- unique_input$MAF_Chinese[idx_m]
          test_rs[j, "MAF_Indian"]  <- unique_input$MAF_Indian[idx_m]
          test_rs[j, "MAF_Malay"]   <- unique_input$MAF_Malay[idx_m]
        }

        if (nrow(unique_input) > MIN_SNPS_FOR_COLOC) {
          is_case_control <- sample_size_gwas$Case_control[
            which(sample_size_gwas$Study == trait)]

          if (is_case_control == "Yes") {
            # Case-control GWAS: type "cc" with the case proportion s
            result <- coloc.abf(
              dataset1 = list(pvalues = unique_input$Pvalue_gwas, type = "cc",
                              beta = unique_input$beta_gwas,
                              varbeta = unique_input$varbeta,
                              s = as.numeric(sample_size_gwas$Case_proportion[
                                sample_size_gwas$Study == trait]),
                              snp = unique_input$variant_id),
              dataset2 = list(pvalues = unique_input$p.value_qtl, type = "quant",
                              N = unique(snp_1mb$N_qtl),
                              snp = unique_input$variant_id),
              MAF = unique_input$MAF_Chinese)
          } else if (is_case_control == "No") {
            # Quantitative GWAS
            result <- coloc.abf(
              dataset1 = list(pvalues = unique_input$Pvalue_gwas, type = "quant",
                              beta = unique_input$beta_gwas,
                              varbeta = unique_input$varbeta,
                              N = as.numeric(sample_size_gwas$N[
                                sample_size_gwas$Study == trait]),
                              snp = unique_input$variant_id),
              dataset2 = list(pvalues = unique_input$p.value_qtl, type = "quant",
                              N = unique(snp_1mb$N_qtl),
                              snp = unique_input$variant_id),
              MAF = unique_input$MAF_Chinese)
          }

          test_rs[j, 3:8] <- t(as.data.frame(result$summary))[1, 1:6]

          # Keep the input summary statistics of colocalizing loci for plotting
          if (!is.null(result$summary["PP.H4.abf"]) &&
              as.numeric(result$summary["PP.H4.abf"]) > PP4_SAVE_THRESHOLD) {
            coloc_sumstat_dir <- paste0(output_coloc_file_prefix, ".", trait_short)
            dir.create(coloc_sumstat_dir, recursive = TRUE, showWarnings = FALSE)
            fwrite(unique_input,
                   file.path(coloc_sumstat_dir,
                             paste0(phenotype, ".sum_stat.txt")), sep = "\t")
          }
        }
      }

      test_rs$trait <- trait
      write.table(test_rs,
                  paste0(output_coloc_file_prefix, ".", trait_short, ".tsv"),
                  sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
      message("  - Coloc results written to file.")
    }
    message(">>> Finished trait: ", trait, "\n")
  }
}

# ====== Entry point ======
# Two ways to run this script:
#   * with -i/-o/-c/-d/-N : one explicit input file
#   * with no options     : every cell type and depot of the reference lists,
#     with all paths following the repository layout (-t restricts it to one
#     cell type)
# The GWAS summary statistics are re-read for every cell type.

if (!is.null(args$input_file)) {
  if (is.null(args$output_inter_file) ||
      is.null(args$output_coloc_file_prefix) ||
      is.null(args$qtl_sum_stats_dir) || is.null(args$N)) {
    cat("Error: -i also requires -o, -c, -d and -N. Use -h for usage.\n")
    q(status = 1)
  }
  run_cell_type_depot(args$input_file, args$output_inter_file,
                      args$output_coloc_file_prefix, args$qtl_sum_stats_dir,
                      args$N)
} else {
  model_dir <- "../../data/05_sQTL_mapping"
  coloc_dir <- "../../data/07_colocalization"
  ref_dir   <- "../../reference"
  n_donors  <- 296
  pc_splicing <- 5   # splicing PCs, part of the result file names

  # Genotype PCs used per depot (matches the mapping stage)
  g_for_depot <- function(depot) if (depot == "SAT") 10 else 6

  for (depot in c("SAT", "VAT")) {
    g <- g_for_depot(depot)

    cell_types <- if (!is.null(args$cell_type)) {
      args$cell_type
    } else {
      f <- file.path(ref_dir, paste0("celltypes.QTL.", depot, ".txt"))
      if (!file.exists(f)) stop("Cell-type list not found: ", f)
      read.table(f, sep = "\t", stringsAsFactors = FALSE)[[1]]
    }

    for (ct in cell_types) {
      loci_file <- file.path(
        model_dir, "results",
        sprintf("%s.%s.p%d_g%d.result.lead_variants.fdr", ct, depot, pc_splicing, g))

      if (!file.exists(loci_file) || file.size(loci_file) == 0) {
        message("skip ", ct, " (", depot, "): no loci file")
        next
      }

      message("=== ", ct, " (", depot, ") ===")
      out_dir <- file.path(coloc_dir, "sQTL_and_GWAS", ct, depot)
      qtl_dir <- file.path(model_dir, "cisQTL", ct, depot,
                           sprintf("p%d_g%d", pc_splicing, g))

      run_cell_type_depot(loci_file,
                          file.path(out_dir, "intersect.loci"),
                          file.path(out_dir, "coloc"),
                          qtl_dir,
                          n_donors)
    }
  }
}
