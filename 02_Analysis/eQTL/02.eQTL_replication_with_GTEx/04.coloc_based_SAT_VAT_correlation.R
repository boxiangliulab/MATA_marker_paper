#!/usr/bin/env Rscript
# 04.coloc_based_SAT_VAT_correlation.R
#
# Purpose:
#   Use coloc-defined categories (Shared/SAT_specific/VAT_specific) for Adipocyte,
#   then for each gene pull ALL SNPs from SAT and VAT cis files, take the
#   intersection of tested SNPs, and compute Pearson r between BETA_SAT and
#   BETA_VAT across all SNP×gene pairs within each category.
#
# Input:
#   Coloc results: /home/users/nus/e1101920/.../eQTL_SAT_VAT_coloc/Adipocyte/
#   Cis files:     .../03.formal_calling/ancestry_specific_version/{SAT,VAT}/Adipocyte/
#
# Output (OUTPUT_DIR):
#   coloc_category_all_snp_correlation.tsv   -- per-category Pearson r
#   coloc_category_per_gene_r.tsv            -- per-gene r (distribution)
#   coloc_category_correlation_boxplot.pdf   -- boxplot of per-gene r
#   coloc_category_correlation_scatter.pdf   -- scatter sample per category

library(data.table)
library(ggplot2)

## =========================================================
## Paths
## =========================================================
COLOC_DIR <- "/home/users/nus/e1101920/project/adipose/05_model/depot_eQTL/results/eQTL_SAT_VAT_coloc/Adipocyte"
CIS_ROOT  <- "/data/projects/11003054/e1101943/project/01.adipose_sqtl/data/09.eQTL/03.Saige_QTL/02.output/03.formal_calling/ancestry_specific_version"
OUTPUT_DIR <- "/data/projects/11003054/e1101943/project/01.adipose_sqtl/data/09.eQTL/03.Saige_QTL/02.output/04.replication/coloc_based_correlation"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

## =========================================================
## Read coloc gene lists
## =========================================================
shared_genes      <- fread(file.path(COLOC_DIR, "colocalized.tsv"), header = TRUE)$gene
sat_specific_genes <- fread(file.path(COLOC_DIR, "SAT_specific.tsv"), header = TRUE)$gene
vat_specific_genes <- fread(file.path(COLOC_DIR, "VAT_specific.tsv"), header = TRUE)$gene

message("Shared (colocalized): ", length(shared_genes))
message("SAT_specific: ",         length(sat_specific_genes))
message("VAT_specific: ",         length(vat_specific_genes))

gene_list <- list(
  Shared       = shared_genes,
  SAT_specific = sat_specific_genes,
  VAT_specific = vat_specific_genes
)

## =========================================================
## Helper: read cis file, return data.table(MarkerID, BETA)
## =========================================================
read_cis <- function(cis_root, depot, gene) {
  f <- file.path(cis_root, depot, "Adipocyte", paste0(gene, "_cis"))
  if (!file.exists(f)) return(NULL)
  fread(f, select = c("MarkerID", "BETA"), showProgress = FALSE)
}

## =========================================================
## Main: for each category, collect all intersection SNPs
## =========================================================
per_gene_list   <- list()
pooled_list     <- list()

for (cat_name in names(gene_list)) {
  genes <- gene_list[[cat_name]]
  message("\n=== ", cat_name, " (", length(genes), " genes) ===")

  cat_rows <- rbindlist(lapply(seq_along(genes), function(i) {
    g <- genes[i]
    if (i %% 200 == 0) message("  ", i, "/", length(genes))

    sat_dt <- read_cis(CIS_ROOT, "SAT", g)
    vat_dt <- read_cis(CIS_ROOT, "VAT", g)

    if (is.null(sat_dt) || is.null(vat_dt)) return(NULL)
    if (nrow(sat_dt) == 0 || nrow(vat_dt) == 0) return(NULL)

    # Intersection of tested SNPs
    merged <- merge(
      sat_dt[, .(MarkerID, BETA_SAT = BETA)],
      vat_dt[, .(MarkerID, BETA_VAT = BETA)],
      by = "MarkerID"
    )
    if (nrow(merged) == 0) return(NULL)

    merged[, gene := g]
    merged
  }), fill = TRUE)

  if (is.null(cat_rows) || nrow(cat_rows) == 0) next

  cat_rows[, category := cat_name]
  pooled_list[[cat_name]] <- cat_rows

  # Per-gene r
  gene_r <- cat_rows[, {
    if (.N >= 3) {
      ct <- cor.test(BETA_SAT, BETA_VAT)
      list(
        pearson_r      = unname(ct$estimate),
        pearson_p      = ct$p.value,
        n_snps         = .N,
        pct_concordant = round(100 * sum(sign(BETA_SAT) == sign(BETA_VAT)) / .N, 1)
      )
    } else {
      list(pearson_r = NA_real_, pearson_p = NA_real_,
           n_snps = .N, pct_concordant = NA_real_)
    }
  }, by = gene]
  gene_r[, category := cat_name]
  per_gene_list[[cat_name]] <- gene_r

  message("  Total intersection SNPs: ", nrow(cat_rows))
  message("  Genes with data: ", nrow(gene_r))
}

## =========================================================
## Pooled correlation per category
## =========================================================
all_pooled <- rbindlist(pooled_list, fill = TRUE)

pool_stat <- all_pooled[, {
  ct <- cor.test(BETA_SAT, BETA_VAT)
  list(
    pearson_r      = unname(ct$estimate),
    pearson_p      = ct$p.value,
    n_snp_pairs    = .N,
    n_genes        = uniqueN(gene),
    pct_concordant = round(100 * sum(sign(BETA_SAT) == sign(BETA_VAT)) / .N, 1)
  )
}, by = category]

pool_stat[, category := factor(category,
                                levels = c("Shared", "SAT_specific", "VAT_specific"))]
setorder(pool_stat, category)

fwrite(pool_stat,
       file = file.path(OUTPUT_DIR, "coloc_category_all_snp_correlation.tsv"),
       sep = "\t")
message("\nPooled correlation table:")
print(pool_stat)

## =========================================================
## Per-gene r table
## =========================================================
per_gene_dt <- rbindlist(per_gene_list, fill = TRUE)
per_gene_dt[, category := factor(category,
                                   levels = c("Shared", "SAT_specific", "VAT_specific"))]
fwrite(per_gene_dt,
       file = file.path(OUTPUT_DIR, "coloc_category_per_gene_r.tsv"),
       sep = "\t")

## =========================================================
## Boxplot of per-gene r
## =========================================================
p_box <- ggplot(per_gene_dt[!is.na(pearson_r)],
                aes(x = category, y = pearson_r, fill = category)) +
  geom_boxplot(outlier.size = 0.5, width = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  scale_fill_manual(values = c(
    "Shared"       = "black",
    "SAT_specific" = "#E69F00",
    "VAT_specific" = "#56B4E9"
  )) +
  labs(
    x     = NULL,
    y     = "Per-gene Pearson r (BETA_SAT vs BETA_VAT, all cis-SNPs)",
    title = "SAT-VAT eQTL effect-size concordance by coloc category\n(Adipocyte, all cis-SNPs)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none"
  )

ggsave(
  filename = file.path(OUTPUT_DIR, "coloc_category_correlation_boxplot.pdf"),
  plot = p_box, width = 5, height = 5, useDingbats = FALSE
)

## =========================================================
## Scatter plot (sample 500 SNPs per category for visualisation)
## =========================================================
set.seed(42)
scatter_dt <- all_pooled[, .SD[sample(.N, min(.N, 500))], by = category]

p_scatter <- ggplot(scatter_dt,
                    aes(x = BETA_SAT, y = BETA_VAT, color = category)) +
  geom_point(alpha = 0.3, size = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  facet_wrap(~category, nrow = 1) +
  scale_color_manual(values = c(
    "Shared"       = "black",
    "SAT_specific" = "#E69F00",
    "VAT_specific" = "#56B4E9"
  )) +
  labs(
    x     = "BETA_SAT",
    y     = "BETA_VAT",
    title = "SAT vs VAT effect sizes by coloc category (sampled SNPs)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none"
  )

ggsave(
  filename = file.path(OUTPUT_DIR, "coloc_category_correlation_scatter.pdf"),
  plot = p_scatter, width = 10, height = 4, useDingbats = FALSE
)

message("\nDone. Output in: ", OUTPUT_DIR)
