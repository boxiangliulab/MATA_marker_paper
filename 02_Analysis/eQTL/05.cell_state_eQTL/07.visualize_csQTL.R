# This script is to visualize the cell state QTLs
# Yihan Tong
# 2026.3.11

# Project root is 8 levels up from this script
args0 <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args0[grep("--file=", args0)])
project_root <- normalizePath(file.path(dirname(normalizePath(script_path)),
                                        paste(rep("..", 8), collapse = "/")))

base_data <- file.path(project_root, "data/09.eQTL/03.Saige_QTL/02.output/11.cell_state_eQTL/02.adipogenesis")

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(ggplot2)
  library(patchwork)
})
share_csqtl <- read.table(file.path(base_data, "SAT_pc2_VAT_pc3.sig_cseQTL.shared.FDR0.05.tsv"), sep = "\t", header= T)
genotype <- fread("/home/project/11003054/e1353486/project/adipose/fullgenotype/Ancestry_specific/00.merge_genotype/genotype_adipose_after_imputation_296_merged_genome_rename_rsID_neur.tsv",sep = "\t")
share_csqtl_clean <- share_csqtl %>%
  rename(
    snp_lead_sat = snp_lead.x,
    snp_lead_vat = snp_lead.y,
    term_sat     = term.x,
    term_vat     = term.y,
    estimate_sat = estimate.x,
    estimate_vat = estimate.y,
    se_sat       = se.x,
    se_vat       = se.y,
    pvalue_sat   = pvalue.x,
    pvalue_vat   = pvalue.y,
    pPC_sat      = pPC.x,
    pPC_vat      = pPC.y,
    depot_sat    = depot.x,
    depot_vat    = depot.y,
    FDR_sat      = FDR.x,
    FDR_vat      = FDR.y,
    eQTL_sat     = eQTL.x,
    eQTL_vat     = eQTL.y
  ) %>%
  filter(snp_id == snp_lead_sat)
share_csqtl_clean <- as.data.table(share_csqtl_clean)
# read phenotype
read_data <- function(depot, pc_quantile) {
  pheno_dir <- file.path(base_data, depot, "01.phenotype_for_separate_run",
                         paste0("phenotype_", pc_quantile, ".tsv"))
  data <- fread(pheno_dir, quote="\t", header=T)
}

depot <- "SAT"   # "SAT" or "VAT"
# Three quantiles to iterate over
pc_quantiles <- c("high", "low", "median")

# Output directory
outdir <- file.path(base_data, "shared_cs_eQTL_plot")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Root directory for cis result files
cis_root <- file.path(base_data, depot, "output", "separate_test_result")

# Genotype colors
geno_cols <- c("0" = "#BFDFD2", "1" = "#DBCB92", "2" = "#ED8D5A")

# Number of metadata columns
meta_n <- if (toupper(depot) == "SAT") 18 else 16

# FDR column in share_csqtl_clean
fdr_col <- paste0("FDR_", tolower(depot))

# =========================================================
# Read phenotype data for one quantile
# =========================================================
read_data <- function(depot, pc_quantile) {
  pheno_file <- file.path(base_data, depot, "01.phenotype_for_separate_run",
                          paste0("phenotype_", pc_quantile, ".tsv"))
  dt <- fread(pheno_file, header = TRUE)
  return(dt)
}

# =========================================================
# Find cis file for one gene within one quantile
# =========================================================
find_cis_file <- function(cis_root, pc_quantile, gene) {
  dir_path <- file.path(cis_root, pc_quantile)
  if (!dir.exists(dir_path)) return(NA_character_)

  files <- list.files(dir_path, full.names = TRUE)
  if (length(files) == 0) return(NA_character_)

  bns <- basename(files)

  patterns_exact <- c(
    paste0("^", gene, "_cis(\\..+)?$")
  )

  for (pat in patterns_exact) {
    hit <- files[grepl(pat, bns)]
    if (length(hit) >= 1) return(hit[1])
  }

  return(NA_character_)
}
# =========================================================
# Extract BETA for a gene-SNP pair from cis result
# =========================================================
read_beta_from_cis <- function(cis_root, pc_quantile, gene, snp_id) {
  cis_file <- find_cis_file(cis_root, pc_quantile, gene)
  if (is.na(cis_file) || !file.exists(cis_file)) return(NA_real_)

  dt <- tryCatch(
    fread(cis_file),
    error = function(e) NULL
  )
  if (is.null(dt) || nrow(dt) == 0) return(NA_real_)

  marker_col <- c("MarkerID", "SNP", "snp_id")[c("MarkerID", "SNP", "snp_id") %in% colnames(dt)][1]
  beta_col   <- c("BETA", "Beta", "beta")[c("BETA", "Beta", "beta") %in% colnames(dt)][1]

  if (is.na(marker_col) || is.na(beta_col)) return(NA_real_)

  hit <- dt[get(marker_col) == snp_id]
  if (nrow(hit) == 0) return(NA_real_)

  return(as.numeric(hit[[beta_col]][1]))
}

# =========================================================
# Plot one gene-SNP-quantile combination
# =========================================================
make_one_plot <- function(dt_sub, depot) {
  dt_sub <- copy(dt_sub)
  dt_plot <- dt_sub[expr > 0]

  if (nrow(dt_plot) == 0) return(NULL)

  gene_now <- unique(dt_sub$gene_id)[1]
  snp_now  <- unique(dt_sub$snp_id)[1]
  q_now    <- unique(dt_sub$pc_quantile)[1]
  fdr_now  <- unique(dt_sub$FDR)[1]
  beta_now <- unique(dt_sub$BETA)[1]

  y_max <- max(dt_plot$expr, na.rm = TRUE)
  y_min <- min(dt_plot$expr, na.rm = TRUE)
  y_pos <- y_max + (y_max - y_min) * 0.12
  if (!is.finite(y_pos)) y_pos <- y_max + 0.5

  label_txt <- paste0(
    "BETA = ", ifelse(is.na(beta_now), "NA", sprintf("%.4f", beta_now)),
    "\n",
    "FDR_", tolower(depot), " = ", ifelse(is.na(fdr_now), "NA", sprintf("%.3g", fdr_now)),
    "\n",
    "pc_quantile = ", q_now
  )

  p <- ggplot(dt_plot, aes(x = genotype, y = expr, fill = genotype)) +
    geom_violin(trim = TRUE, color = "black", linewidth = 0.35, alpha = 0.9) +
    geom_boxplot(width = 0.16, outlier.shape = NA, fill = "white", color = "black", linewidth = 0.35) +
    geom_jitter(width = 0.08, size = 1.2, alpha = 0.8, color = "black") +
    scale_fill_manual(values = geno_cols, drop = FALSE) +
    labs(
      title = paste0(gene_now, " | ", snp_now),
      subtitle = paste0(depot, " | pc_quantile: ", q_now, " | zero-expression donors removed"),
      x = "Genotype",
      y = "Pseudobulk expression: log2(CPM + 1)"
    ) +
    annotate("text", x = 2.8, y = y_pos, label = label_txt, hjust = 1, vjust = 1, size = 4) +
    coord_cartesian(ylim = c(y_min, y_pos + (y_max - y_min) * 0.05), clip = "off") +
    theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 11),
      axis.text.x = element_text(color = "black"),
      panel.grid = element_blank(),
      legend.position = "none"
    )

  return(p)
}

# =========================================================
# Process each pc_quantile
# =========================================================
all_plots <- list()

for (pc_q in pc_quantiles) {

  message("Processing pc_quantile: ", pc_q)

  # -------------------------------------------------------
  # 1. Read phenotype file
  # -------------------------------------------------------
  data <- read_data(depot = depot, pc_quantile = pc_q)
  data <- as.data.table(data)

  if (!"donor_id" %in% colnames(data)) {
    stop(paste0("donor_id column is missing in phenotype_", pc_q, ".tsv"))
  }
  if (ncol(data) <= meta_n) {
    stop(paste0("Not enough columns in phenotype_", pc_q, ".tsv for meta/gene split."))
  }

  # -------------------------------------------------------
  # 2. Split metadata and gene count matrix
  # -------------------------------------------------------
  meta_dt <- data[, 1:meta_n, with = FALSE]
  gene_dt <- data[, (meta_n + 1):ncol(data), with = FALSE]

  for (j in seq_along(gene_dt)) {
    set(gene_dt, j = j, value = as.numeric(gene_dt[[j]]))
  }

  genes <- colnames(gene_dt)

  # -------------------------------------------------------
  # 3. Pseudobulk by donor_id
  # -------------------------------------------------------
  pb_counts_dt <- cbind(meta_dt[, .(donor_id)], gene_dt)
  pb_counts_dt <- pb_counts_dt[, lapply(.SD, sum, na.rm = TRUE), by = donor_id, .SDcols = genes]

  pb_counts_mat <- as.matrix(pb_counts_dt[, ..genes])
  rownames(pb_counts_mat) <- pb_counts_dt$donor_id
  storage.mode(pb_counts_mat) <- "numeric"

  fwrite(
    data.table(donor_id = rownames(pb_counts_mat), as.data.table(pb_counts_mat)),
    file.path(outdir, paste0(depot, "_pcq", pc_q, "_donor_by_gene_pseudobulk_counts.tsv.gz")),
    sep = "\t"
  )

  # -------------------------------------------------------
  # 4. Compute log2(CPM + 1) using edgeR
  # -------------------------------------------------------
  dge <- DGEList(counts = t(pb_counts_mat))   # genes x donors
  cpm_mat <- cpm(dge, log = FALSE, normalized.lib.sizes = TRUE)
  log2cpm_mat <- log2(cpm_mat + 1)

  log2cpm_donor_gene <- t(log2cpm_mat)

  fwrite(
    data.table(donor_id = rownames(log2cpm_donor_gene), as.data.table(log2cpm_donor_gene)),
    file.path(outdir, paste0(depot, "_pcq", pc_q, "_donor_by_gene_log2CPMplus1.tsv.gz")),
    sep = "\t"
  )

  # -------------------------------------------------------
  # 5. Keep gene-SNP pairs for this quantile from share_csqtl_clean
  # -------------------------------------------------------
 plot_df0 <- unique(
  share_csqtl_clean[
    !is.na(gene_id) &
      !is.na(snp_id) &
      gene_id %in% colnames(log2cpm_donor_gene),
    .(
      gene_id = as.character(gene_id),
      snp_id = as.character(snp_id),
      pc_quantile = as.character(pc_q),
      FDR = as.numeric(get(fdr_col))
    )
  ]
)
  if (nrow(plot_df0) == 0) {
    message("No genes found in share_csqtl_clean for pc_quantile ", pc_q)
    next
  }

  # -------------------------------------------------------
  # 6. Convert genotype matrix to long format
  # -------------------------------------------------------
  geno_long <- melt(
    genotype,
    id.vars = "SNP",
    variable.name = "donor_id",
    value.name = "genotype"
  )
  geno_long[, donor_id := as.character(donor_id)]
  geno_long[, SNP := as.character(SNP)]
  geno_long[, genotype := as.character(genotype)]

  # -------------------------------------------------------
  # 7. Convert expression matrix to long format
  # -------------------------------------------------------
  expr_long <- as.data.table(log2cpm_donor_gene, keep.rownames = "donor_id")
  expr_long <- melt(
    expr_long,
    id.vars = "donor_id",
    variable.name = "gene_id",
    value.name = "expr"
  )
  expr_long[, donor_id := as.character(donor_id)]
  expr_long[, gene_id := as.character(gene_id)]

  # -------------------------------------------------------
  # 8. Merge gene-SNP list, expression, and genotype
  # -------------------------------------------------------
  plot_long <- merge(
    plot_df0,
    expr_long,
    by = "gene_id",
    allow.cartesian = TRUE
  )

  plot_long <- merge(
    plot_long,
    geno_long,
    by.x = c("snp_id", "donor_id"),
    by.y = c("SNP", "donor_id"),
    all.x = TRUE
  )

  plot_long <- plot_long[!is.na(genotype)]
  plot_long[, genotype := factor(genotype, levels = c("0", "1", "2"))]

  # -------------------------------------------------------
  # 9. Read BETA from cis result
  # -------------------------------------------------------
  beta_info <- unique(plot_long[, .(gene_id, snp_id, pc_quantile, FDR)])

  beta_info[, BETA := mapply(
    FUN = read_beta_from_cis,
    pc_quantile = pc_quantile,
    gene = gene_id,
    snp_id = snp_id,
    MoreArgs = list(cis_root = cis_root)
  )]

  plot_long <- merge(
    plot_long,
    beta_info,
    by = c("gene_id", "snp_id", "pc_quantile", "FDR"),
    all.x = TRUE
  )
  plot_long$depot=depot
  fwrite(
    plot_long,
    file.path(outdir, paste0(depot, "_pcq", pc_q, "_gene_snp_donor_expression_genotype_merged.tsv.gz")),
    sep = "\t"
  )

  # -------------------------------------------------------
  # 10. Generate plots
  # -------------------------------------------------------
  plot_keys <- unique(plot_long[, .(gene_id, snp_id, pc_quantile)])

  for (i in seq_len(nrow(plot_keys))) {
    g <- plot_keys$gene_id[i]
    s <- plot_keys$snp_id[i]
    q <- plot_keys$pc_quantile[i]

    dt_sub <- plot_long[gene_id == g & snp_id == s & pc_quantile == q]
    if (nrow(dt_sub) == 0) next

    p <- make_one_plot(dt_sub, depot = depot)

    plot_name <- paste(g, s, q, sep = "__")
    all_plots[[plot_name]] <- p

    ggsave(
      file.path(outdir,q, paste0(depot, "_", g, "__", s, "__pcq_", q, ".pdf")),
      p,
      width = 4,
      height = 5
    )
  }
}

# =========================================================
# Save combined plot
# =========================================================
if (length(all_plots) > 0) {
  combined <- wrap_plots(all_plots, ncol = 3)

  ggsave(
    file.path(outdir, paste0(depot, "_all_genotype_violin_plots.pdf")),
    combined,
    width = 15.6,
    height = 5.2 * ceiling(length(all_plots) / 3)
  )
}

message("All done.")
message("Output directory: ", outdir)




