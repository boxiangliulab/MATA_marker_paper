# Per-gene pseudobulk expression boxplot ~ BMI group
# for nominally significant genes (lm_p < 0.05) from 03b heatmap analysis
# Only donors with >10 committed ASPC cells (already filtered in 03a pseudobulk)
#
# Input:
#   data/30.adipogenesis/03.committed_aspc_pseudobulk/{SAT,VAT}.pseudobulk_counts.tsv
#   data/30.adipogenesis/03.committed_aspc_gene_BMI_stats.tsv
# Output:
#   data/30.adipogenesis/05.committed_aspc_fibrosis_expression.pdf
#   data/30.adipogenesis/05.committed_aspc_adipogenesis_expression.pdf

suppressPackageStartupMessages({
  library(edgeR)
  library(msigdbr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(data.table)
  library(cowplot)
})

# Project root is 7 levels up from this script
args0 <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args0[grep("--file=", args0)])
project_root <- normalizePath(file.path(dirname(normalizePath(script_path)),
                                        paste(rep("..", 7), collapse = "/")))

IN_DIR     <- file.path(project_root, "data/30.adipogenesis/03.committed_aspc_pseudobulk")
STATS_FILE <- file.path(project_root, "data/30.adipogenesis/03.committed_aspc_gene_BMI_stats.tsv")
META_FILE  <- file.path(project_root, "data/08.integration_snRNA_object/01.freeze_version/metadata_after_libraryQC_after_ctQC_freeze.tsv")
OUT_DIR    <- file.path(project_root, "data/30.adipogenesis")

DEPOT_COLORS <- c(SAT = "#67B3AD", VAT = "#F3AF44")
BMI_LEVELS   <- c("Non-obese", "Obese", "Severely obese", "Morbidly obese")
ASIAN_RACES  <- c("Chinese", "Malay", "Indian")

# ── Gene sets ─────────────────────────────────────────────────────────────────
msig_all <- msigdbr(species = "Homo sapiens")
get_msig <- function(gs) unique(msig_all$gene_symbol[msig_all$gs_name == gs])

GENESETS <- list(
  eighteen_pro_fibrotic_genes       = c("COL1A1","COL3A1","COL6A1","COL6A3","LUM","TNC",
                         "COL5A1","LOX","FN1","LOXL1","LOXL2","CCN2",
                         "TGFB1","VCAN","LRP1","THBS1","THY1","ACTA2"),
  NABA_Collagens     = get_msig("NABA_COLLAGENS"),
  Collagen_Formation = get_msig("REACTOME_COLLAGEN_FORMATION"),
  HALLMARK           = get_msig("HALLMARK_ADIPOGENESIS"),
  REACTOME           = get_msig("REACTOME_ADIPOGENESIS"),
  GOBP_WFC           = get_msig("GOBP_WHITE_FAT_CELL_DIFFERENTIATION"),
  KEGG_PPAR          = get_msig("KEGG_PPAR_SIGNALING_PATHWAY")
)

fibrosis_genes   <- unique(unlist(GENESETS[c("eighteen_pro_fibrotic_genes","NABA_Collagens","Collagen_Formation")]))
adipogenesis_genes <- unique(unlist(GENESETS[c("HALLMARK","REACTOME","GOBP_WFC","KEGG_PPAR")]))

# ── Load pseudobulk counts → log1p(CPM) ───────────────────────────────────────
load_logcpm <- function(depot) {
  counts <- read.table(file.path(IN_DIR, paste0(depot, ".pseudobulk_counts.tsv")),
                       sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
  dge    <- calcNormFactors(DGEList(counts = as.matrix(counts)))
  cpm_mat <- cpm(dge, normalized.lib.sizes = TRUE, log = FALSE)
  log1p(cpm_mat)
}

logcpm_list        <- lapply(c("SAT","VAT"), load_logcpm)
names(logcpm_list) <- c("SAT","VAT")

# ── Metadata ──────────────────────────────────────────────────────────────────
meta <- fread(META_FILE, data.table = FALSE, na.strings = c("", "NA")) %>%
  filter(Race %in% ASIAN_RACES, !is.na(BMI)) %>%
  mutate(
    BMI       = as.numeric(BMI),
    BMI_group = case_when(
      BMI <  27.5 ~ "Non-obese",
      BMI <  32.5 ~ "Obese",
      BMI <  37.5 ~ "Severely obese",
      TRUE        ~ "Morbidly obese"
    ),
    BMI_group     = factor(BMI_group, levels = BMI_LEVELS),
    BMI_group_num = as.numeric(factor(BMI_group, levels = BMI_LEVELS))
  ) %>%
  select(donor_id, BMI_group, BMI_group_num) %>%
  distinct(donor_id, .keep_all = TRUE)

# ── Load stats, filter nominal p < 0.05 ───────────────────────────────────────
stats <- fread(STATS_FILE, data.table = FALSE)

sig_genes <- stats %>%
  filter(lm_p < 0.05) %>%
  pull(gene) %>%
  unique()

cat(sprintf("Nominally significant genes (lm_p<0.05 in SAT or VAT): %d\n", length(sig_genes)))

# ── Build long expression data frame ──────────────────────────────────────────
expr_long <- bind_rows(lapply(c("SAT","VAT"), function(dp) {
  mat  <- logcpm_list[[dp]]
  keep <- intersect(sig_genes, rownames(mat))
  mat[keep, , drop = FALSE] %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "donor_id", values_to = "log1p_cpm") %>%
    mutate(Depot = dp)
})) %>%
  inner_join(meta, by = "donor_id") %>%
  mutate(Depot = factor(Depot, levels = c("SAT","VAT")))

# ── Stats annotation per gene × depot ─────────────────────────────────────────
star_label <- function(p) case_when(
  p < 0.001 ~ "***",
  p < 0.01  ~ "**",
  p < 0.05  ~ "*",
  p < 0.1   ~ ".",
  TRUE      ~ ""
)

# one combined label per gene: SAT line + VAT line
stats_annot <- stats %>%
  filter(gene %in% sig_genes) %>%
  mutate(stars = star_label(lm_p),
         line  = paste0(toupper(depot), ": beta=",
                        formatC(lm_beta, format = "e", digits = 1),
                        " p=", formatC(lm_p, format = "e", digits = 1),
                        " ", stars)) %>%
  group_by(gene) %>%
  summarise(label = paste(line, collapse = "\n"), .groups = "drop")

# ── Plot function: one page of genes ──────────────────────────────────────────
plot_gene_grid <- function(genes, title, outfile, ncol_plot = 4) {
  df   <- expr_long %>% filter(gene %in% genes)
  ann  <- stats_annot %>% filter(gene %in% genes)

  # y position for annotation: per gene (use global max across both depots)
  ymax_df <- df %>%
    group_by(gene) %>%
    summarise(ymax = max(log1p_cpm, na.rm = TRUE), .groups = "drop")
  ann <- ann %>% left_join(ymax_df, by = "gene")

  plts <- lapply(genes, function(g) {
    gdf  <- df %>% filter(gene == g)
    adf  <- ann %>% filter(gene == g)
    ggplot(gdf, aes(x = BMI_group_num, y = log1p_cpm, color = Depot, fill = Depot)) +
      geom_boxplot(aes(group = interaction(BMI_group_num, Depot)),
                   outlier.shape = NA, alpha = 0.25, linewidth = 0.35,
                   width = 0.33, position = position_dodge(0.6)) +
      geom_jitter(size = 0.5, alpha = 0.5,
                  position = position_jitterdodge(jitter.width = 0.12,
                                                  dodge.width  = 0.6)) +
      geom_smooth(method = "lm", se = TRUE, alpha = 0.15,
                  linewidth = 0.7, fullrange = FALSE) +
      annotate("text", x = 1, y = adf$ymax * 1.02,
               label = adf$label, size = 1.8, hjust = 0, vjust = 1,
               color = "grey20", lineheight = 1.1) +
      scale_fill_manual(values  = DEPOT_COLORS) +
      scale_color_manual(values = DEPOT_COLORS) +
      scale_x_continuous(breaks = 1:4, labels = BMI_LEVELS) +
      coord_fixed(ratio = 3 / diff(range(gdf$log1p_cpm, na.rm = TRUE))) +
      labs(x = NULL, y = "log1p(CPM)", title = g) +
      theme_classic(base_size = 8) +
      theme(
        axis.text.x     = element_text(angle = 30, hjust = 1, size = 7),
        plot.title      = element_text(face = "bold.italic", size = 8),
        legend.position = "none"
      )
  })

  n_pages <- ceiling(length(genes) / (ncol_plot * 4))
  per_page <- ncol_plot * 4

  pdf(outfile, width = ncol_plot * 2.4, height = ceiling(per_page / ncol_plot) * 2.4)
  for (i in seq(1, length(genes), per_page)) {
    chunk <- plts[i:min(i + per_page - 1, length(genes))]
    print(plot_grid(plotlist = chunk, ncol = ncol_plot))
  }
  # legend page
  leg_p <- ggplot(df %>% filter(gene == genes[1]),
                  aes(x = BMI_group, y = log1p_cpm, fill = Depot, color = Depot)) +
    geom_boxplot() +
    scale_fill_manual(values  = DEPOT_COLORS) +
    scale_color_manual(values = DEPOT_COLORS) +
    theme_classic(base_size = 9) +
    theme(legend.position = "bottom")
  print(get_legend(leg_p))
  dev.off()
  message("Done: ", outfile)
}

# ── Fibrosis significant genes ────────────────────────────────────────────────
fib_sig <- intersect(sig_genes, fibrosis_genes)
cat(sprintf("Fibrosis nominal sig genes: %d\n", length(fib_sig)))

# order by SAT lm_beta descending (same as heatmap)
fib_order <- stats %>%
  filter(gene %in% fib_sig, depot == "SAT") %>%
  arrange(desc(lm_beta)) %>% pull(gene) %>%
  union(fib_sig)

if (length(fib_sig) > 0)
  plot_gene_grid(fib_order, "Fibrosis genes",
                 file.path(OUT_DIR, "05.committed_aspc_fibrosis_expression.pdf"))

# ── Adipogenesis significant genes ────────────────────────────────────────────
adipo_sig <- intersect(sig_genes, adipogenesis_genes)
cat(sprintf("Adipogenesis nominal sig genes: %d\n", length(adipo_sig)))

adipo_order <- stats %>%
  filter(gene %in% adipo_sig, depot == "SAT") %>%
  arrange(lm_beta) %>% pull(gene) %>%
  union(adipo_sig)

if (length(adipo_sig) > 0)
  plot_gene_grid(adipo_order, "Adipogenesis genes",
                 file.path(OUT_DIR, "05.committed_aspc_adipogenesis_expression.pdf"))

cat("\nDone.\n")
