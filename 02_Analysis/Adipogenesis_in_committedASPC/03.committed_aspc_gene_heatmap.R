# Per-gene heatmap for collagen/fibrosis and adipogenesis gene sets
# in committed preadipocytes ~ BMI group
#
# Adapted from 29.adipogenesis/06c.preadipocyte_collagen_genes_BMI.R
# Gene sets: same as 02 (no WANG_PPARG)
#
# Input:
#   data/30.adipogenesis/03.committed_aspc_pseudobulk/{SAT,VAT}.pseudobulk_counts.tsv
# Output:
#   data/30.adipogenesis/03.committed_aspc_fibrosis_heatmap.pdf
#   data/30.adipogenesis/03.committed_aspc_adipogenesis_heatmap.pdf
#   data/30.adipogenesis/03.committed_aspc_gene_BMI_stats.tsv

suppressPackageStartupMessages({
  library(edgeR)
  library(msigdbr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(ggpattern)
  library(data.table)
  library(patchwork)
})

IN_DIR    <- "/data/projects/11003054/e1101943/project/01.adipose_sqtl/data/30.adipogenesis/03.committed_aspc_pseudobulk"
META_FILE <- "/data/projects/11003054/e1101943/project/01.adipose_sqtl/data/08.integration_snRNA_object/01.freeze_version/metadata_after_libraryQC_after_ctQC_freeze.tsv"
OUT_DIR   <- "/data/projects/11003054/e1101943/project/01.adipose_sqtl/data/30.adipogenesis"

DEPOT_COLORS <- c(SAT = "#67B3AD", VAT = "#F3AF44")

# ── Gene sets ─────────────────────────────────────────────────────────────────
cat("Fetching MSigDB gene sets...\n")
msig_all <- msigdbr(species = "Homo sapiens")
get_msig <- function(gs) unique(msig_all$gene_symbol[msig_all$gs_name == gs])

GENESETS <- list(
  # Fibrosis
  Eisinger2023       = c("COL1A1","COL3A1","COL6A1","COL6A3","LUM","TNC",
                         "COL5A1","LOX","FN1","LOXL1","LOXL2","CCN2",
                         "TGFB1","VCAN","LRP1","THBS1","THY1","ACTA2"),
  NABA_Collagens     = get_msig("NABA_COLLAGENS"),
  Collagen_Formation = get_msig("REACTOME_COLLAGEN_FORMATION"),
  # Adipogenesis
  HALLMARK           = get_msig("HALLMARK_ADIPOGENESIS"),
  REACTOME           = get_msig("REACTOME_ADIPOGENESIS"),
  GOBP_WFC           = get_msig("GOBP_WHITE_FAT_CELL_DIFFERENTIATION"),
  KEGG_PPAR          = get_msig("KEGG_PPAR_SIGNALING_PATHWAY")
)

GS_GROUPS <- list(
  Fibrosis     = c("Eisinger2023", "NABA_Collagens", "Collagen_Formation"),
  Adipogenesis = c("HALLMARK", "REACTOME", "GOBP_WFC", "KEGG_PPAR")
)

# ── Metadata ──────────────────────────────────────────────────────────────────
meta <- fread(META_FILE, data.table = FALSE, na.strings = c("", "NA")) %>%
  filter(Race %in% c("Chinese", "Malay", "Indian"), !is.na(BMI)) %>%
  mutate(
    BMI           = as.numeric(BMI),
    T2DM          = ifelse(Diabetic == "T2DM", 1, 0),
    BMI_group     = case_when(
      BMI <  27.5 ~ "Non-obese",
      BMI <  32.5 ~ "Obese",
      BMI <  37.5 ~ "Severely obese",
      TRUE        ~ "Morbidly obese"
    ),
    BMI_group     = factor(BMI_group,
                           levels = c("Non-obese","Obese",
                                      "Severely obese","Morbidly obese")),
    BMI_group_num = as.numeric(BMI_group)
  ) %>%
  select(donor_id, BMI, BMI_group, BMI_group_num, Age, Gender, Race, T2DM) %>%
  distinct(donor_id, .keep_all = TRUE)

# ── CPM ───────────────────────────────────────────────────────────────────────
load_cpm <- function(depot) {
  counts <- read.table(file.path(IN_DIR, paste0(depot, ".pseudobulk_counts.tsv")),
                       sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
  dge <- calcNormFactors(DGEList(counts = as.matrix(counts)))
  cpm(dge, normalized.lib.sizes = TRUE)
}

cat("Computing CPM...\n")
cpm_list        <- lapply(c("SAT", "VAT"), load_cpm)
names(cpm_list) <- c("SAT", "VAT")

# ── Helper functions ──────────────────────────────────────────────────────────
fmt_sci     <- function(p) formatC(p, format = "e", digits = 2)
safe_covars <- function(df) {
  covars <- c("Age", "Gender", "T2DM", "Race")
  covars[sapply(covars, function(v) {
    x <- df[[v]]
    if (is.numeric(x)) sd(x, na.rm = TRUE) > 0
    else length(unique(na.omit(x))) >= 2
  })]
}

build_gene_df <- function(genes) {
  bind_rows(lapply(c("SAT", "VAT"), function(depot) {
    mat          <- log1p(cpm_list[[depot]])
    common       <- intersect(colnames(mat), meta$donor_id)
    genes_in     <- intersect(genes, rownames(mat))
    missing      <- setdiff(genes, rownames(mat))
    if (length(missing))
      cat(sprintf("  [%s] not in matrix: %s\n", depot, paste(missing, collapse = ", ")))
    mat[genes_in, common, drop = FALSE] %>%
      as.data.frame() %>%
      rownames_to_column("gene") %>%
      pivot_longer(-gene, names_to = "donor_id", values_to = "expr") %>%
      mutate(depot = depot) %>%
      left_join(meta, by = "donor_id")
  })) %>%
    mutate(depot = factor(depot, levels = c("SAT", "VAT")))
}

run_gene_stats <- function(gene_df) {
  per_depot <- gene_df %>%
    group_by(gene, depot) %>%
    group_modify(~ {
      keep <- safe_covars(.x)
      fmla <- as.formula(paste("expr ~ BMI_group_num +", paste(keep, collapse = " + ")))
      m    <- lm(fmla, data = .x)
      ct   <- summary(m)$coefficients
      sp   <- cor.test(.x$expr, .x$BMI_group_num, method = "spearman", exact = FALSE)
      data.frame(n_donors = nrow(.x),
                 lm_beta  = ct["BMI_group_num", "Estimate"],
                 lm_p     = ct["BMI_group_num", "Pr(>|t|)"],
                 sp_rho   = sp$estimate,
                 sp_p     = sp$p.value)
    }) %>% ungroup()

  interaction <- gene_df %>%
    group_by(gene) %>%
    group_modify(~ {
      if (length(unique(.x$depot)) < 2)
        return(data.frame(int_beta = NA, int_p = NA))
      keep <- safe_covars(.x)
      fmla <- as.formula(paste("expr ~ BMI_group_num * depot +",
                               paste(keep, collapse = " + ")))
      m  <- lm(fmla, data = .x)
      ct <- summary(m)$coefficients
      rw <- "BMI_group_num:depotVAT"
      if (!rw %in% rownames(ct)) return(data.frame(int_beta = NA, int_p = NA))
      data.frame(int_beta = ct[rw, "Estimate"], int_p = ct[rw, "Pr(>|t|)"])
    }) %>% ungroup() %>%
    mutate(int_p_adj = p.adjust(ifelse(is.na(int_p), 1, int_p), method = "BH"))

  per_depot %>%
    left_join(interaction, by = "gene") %>%
    mutate(lm_p_adj = p.adjust(lm_p, method = "BH"))
}

BMI_LEVELS  <- c("Non-obese", "Obese", "Severely obese", "Morbidly obese")
HEAT_COLORS <- c("#313695","#4575b4","#74add1","#abd9e9","#e0f3f8",
                 "#ffffbf","#fee090","#fdae61","#f46d43","#d73027","#a50026")

mean_by_bmi <- function(genes_use) {
  bind_rows(lapply(c("SAT","VAT"), function(depot) {
    mat    <- log1p(cpm_list[[depot]])
    common <- intersect(colnames(mat), meta$donor_id)
    g_in   <- intersect(genes_use, rownames(mat))
    mat[g_in, common, drop = FALSE] %>%
      as.data.frame() %>%
      rownames_to_column("gene") %>%
      pivot_longer(-gene, names_to = "donor_id", values_to = "expr") %>%
      left_join(meta %>% select(donor_id, BMI_group), by = "donor_id") %>%
      mutate(depot = depot)
  })) %>%
    mutate(depot     = factor(depot, levels = c("SAT","VAT")),
           BMI_group = factor(BMI_group, levels = BMI_LEVELS)) %>%
    group_by(gene, depot, BMI_group) %>%
    summarise(mean_logcpm = mean(expr, na.rm = TRUE), .groups = "drop") %>%
    group_by(gene, depot) %>%
    mutate(zscore = as.numeric(scale(mean_logcpm))) %>%
    ungroup()
}

star_label <- function(p) {
  case_when(
    is.na(p)  ~ "",
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ ""
  )
}

# ── draw_heatmap (adapted from 06c) ───────────────────────────────────────────
draw_heatmap <- function(heat_df, sig_genes, annot_df, gene_order,
                         title, out_pdf, gs_labels = NULL,
                         sig_stats = NULL, tile_h = 0.25) {
  gene_fac  <- rev(gene_order)
  gs_orig   <- levels(annot_df$gs_name)
  gs_disp   <- if (!is.null(gs_labels)) {
    v <- gs_labels[gs_orig]; v[is.na(v)] <- gs_orig[is.na(v)]; v
  } else gs_orig
  names(gs_disp) <- gs_orig

  annot_x_levels <- if (!is.null(sig_stats)) c(gs_disp, "Sig.") else gs_disp

  main_df <- heat_df %>%
    filter(gene %in% gene_order) %>%
    mutate(gene  = factor(gene, levels = gene_fac),
           panel = as.character(depot),
           xvar  = factor(as.character(BMI_group), levels = BMI_LEVELS))

  annot_gs_df <- annot_df %>%
    filter(gene %in% gene_order) %>%
    mutate(gene    = factor(gene, levels = gene_fac),
           panel   = "Annot.",
           xvar    = factor(gs_disp[as.character(gs_name)], levels = annot_x_levels),
           pat_val = ifelse(member, "stripe", "none"))

  annot_sig_bg <- if (!is.null(sig_stats))
    tibble(gene  = factor(gene_order, levels = gene_fac),
           panel = "Annot.",
           xvar  = factor("Sig.", levels = annot_x_levels))
  else NULL

  sig_text_df <- if (!is.null(sig_stats)) {
    sig_stats %>%
      filter(gene %in% gene_order) %>%
      mutate(gene  = factor(gene, levels = gene_fac),
             stars = star_label(lm_p),
             panel = "Annot.",
             xvar  = factor("Sig.", levels = annot_x_levels),
             depot = as.character(depot)) %>%
      filter(stars != "")
  } else NULL

  panel_levels <- c("Annot.", "SAT", "VAT")
  main_df$panel     <- factor(main_df$panel,     levels = panel_levels)
  annot_gs_df$panel <- factor(annot_gs_df$panel, levels = panel_levels)

  p <- ggplot() +
    geom_tile(data = main_df,
              aes(x = xvar, y = gene, fill = zscore),
              color = "white", linewidth = 0.25) +
    geom_tile_pattern(
      data            = annot_gs_df,
      aes(x = xvar, y = gene, pattern = pat_val),
      fill            = "white", colour = "black", linewidth = 0.4,
      pattern_colour  = "grey40", pattern_fill = "grey40",
      pattern_angle   = 45, pattern_spacing = 0.05, pattern_density = 0.1,
      inherit.aes     = FALSE
    ) +
    scale_pattern_manual(values = c("stripe" = "stripe", "none" = "none"),
                         guide  = "none")

  if (!is.null(annot_sig_bg)) {
    annot_sig_bg$panel <- factor(annot_sig_bg$panel, levels = panel_levels)
    p <- p + geom_tile(data = annot_sig_bg,
                       aes(x = xvar, y = gene),
                       fill = "white", color = "black", linewidth = 0.4,
                       inherit.aes = FALSE)
  }

  if (!is.null(sig_text_df) && nrow(sig_text_df) > 0) {
    star_size <- max(1.5, 5 * tile_h / 0.3)
    sig_sat   <- sig_text_df %>% filter(depot == "SAT")
    sig_vat   <- sig_text_df %>% filter(depot == "VAT")
    if (nrow(sig_sat) > 0)
      p <- p + geom_text(data = sig_sat,
                         aes(x = xvar, y = gene, label = stars),
                         color = DEPOT_COLORS["SAT"],
                         nudge_y = 0.18, size = star_size, fontface = "bold",
                         inherit.aes = FALSE)
    if (nrow(sig_vat) > 0)
      p <- p + geom_text(data = sig_vat,
                         aes(x = xvar, y = gene, label = stars),
                         color = DEPOT_COLORS["VAT"],
                         nudge_y = -0.18, size = star_size, fontface = "bold",
                         inherit.aes = FALSE)
  }

  p <- p +
    scale_fill_gradientn(colors = HEAT_COLORS, name = "Z-score\n(mean log1p(CPM))",
                         na.value = "white") +
    facet_grid(. ~ panel, scales = "free_x", space = "free_x") +
    labs(x = NULL, y = NULL, title = title) +
    theme_classic(base_size = 9) +
    theme(axis.text.x       = element_text(size = 7, angle = 45, hjust = 1),
          axis.text.y       = element_text(size = 7),
          strip.background  = element_blank(),
          strip.text        = element_text(face = "bold", size = 9),
          legend.key.height = unit(0.8, "cm"),
          plot.title        = element_text(size = 9, face = "bold"),
          panel.spacing     = unit(0.2, "cm"))

  n_genes   <- length(gene_order)
  n_annot_x <- length(annot_x_levels)
  fig_w <- (n_annot_x + 2 * 4) * tile_h + 3.5
  fig_h <- n_genes * tile_h + 2.2
  ggsave(out_pdf, p,
         width = max(5, fig_w), height = max(4, fig_h),
         limitsize = FALSE)
  cat(sprintf("Saved: %s\n", basename(out_pdf)))
}

# ── Run ───────────────────────────────────────────────────────────────────────
all_genes <- unique(unlist(GENESETS))
gene_df   <- build_gene_df(all_genes)
stats     <- run_gene_stats(gene_df)

write.table(stats,
            file.path(OUT_DIR, "03.committed_aspc_gene_BMI_stats.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nFDR < 0.05 hits:\n")
print(stats %>% filter(lm_p_adj < 0.05) %>%
        select(gene, depot, lm_beta, lm_p_adj, sp_rho) %>% arrange(lm_p_adj))

heat_df   <- mean_by_bmi(all_genes)
sig_genes <- stats %>%
  group_by(gene) %>%
  summarise(sig = any(lm_p < 0.05), .groups = "drop") %>%
  filter(sig) %>% pull(gene)
cat(sprintf("\nNominally significant genes (p<0.05 in SAT or VAT): %d\n", length(sig_genes)))

sig_stats <- stats %>% select(gene, depot, lm_p) %>%
  distinct(gene, depot, .keep_all = TRUE)

# ── Fibrosis heatmap ──────────────────────────────────────────────────────────
fib_genes    <- unique(unlist(GENESETS[GS_GROUPS$Fibrosis]))
fib_sig      <- intersect(sig_genes, fib_genes)
fib_order    <- stats %>%
  filter(depot == "SAT", gene %in% fib_sig) %>%
  arrange(desc(lm_beta)) %>% pull(gene) %>%
  union(stats %>% filter(depot == "VAT", gene %in% setdiff(fib_sig, .)) %>%
          arrange(desc(lm_beta)) %>% pull(gene)) %>%
  union(fib_sig)

fib_annot <- tibble(gene = fib_sig) %>%
  mutate(Eisinger2023       = gene %in% GENESETS[["Eisinger2023"]],
         NABA_Collagens     = gene %in% GENESETS[["NABA_Collagens"]],
         Collagen_Formation = gene %in% GENESETS[["Collagen_Formation"]]) %>%
  pivot_longer(-gene, names_to = "gs_name", values_to = "member") %>%
  mutate(gs_name = factor(gs_name, levels = GS_GROUPS$Fibrosis))

fib_gs_labels <- c(
  Eisinger2023       = "Eisinger\n2023",
  NABA_Collagens     = "NABA\nCollagens",
  Collagen_Formation = "Collagen\nFormation"
)

draw_heatmap(heat_df, fib_sig, fib_annot, fib_order,
             "Fibrosis genes ~ BMI\ncommitted preadipocyte (nominal p<0.05)",
             file.path(OUT_DIR, "03.committed_aspc_fibrosis_heatmap.pdf"),
             gs_labels = fib_gs_labels,
             sig_stats = sig_stats %>% filter(gene %in% fib_sig),
             tile_h    = 0.2)

# ── Adipogenesis heatmap ──────────────────────────────────────────────────────
adipo_genes <- unique(unlist(GENESETS[GS_GROUPS$Adipogenesis]))
adipo_sig   <- intersect(sig_genes, adipo_genes)
adipo_order <- stats %>%
  filter(depot == "SAT", gene %in% adipo_sig) %>%
  arrange(lm_beta) %>% pull(gene) %>%
  union(stats %>% filter(depot == "VAT", gene %in% setdiff(adipo_sig, .)) %>%
          arrange(lm_beta) %>% pull(gene)) %>%
  union(adipo_sig)

adipo_annot <- tibble(gene = adipo_sig) %>%
  mutate(HALLMARK        = gene %in% GENESETS[["HALLMARK"]],
         REACTOME        = gene %in% GENESETS[["REACTOME"]],
         GOBP_WFC        = gene %in% GENESETS[["GOBP_WFC"]],
         KEGG_PPAR       = gene %in% GENESETS[["KEGG_PPAR"]]) %>%
  pivot_longer(-gene, names_to = "gs_name", values_to = "member") %>%
  mutate(gs_name = factor(gs_name, levels = GS_GROUPS$Adipogenesis))

adipo_gs_labels <- c(
  HALLMARK        = "HALLMARK",
  REACTOME        = "REACTOME",
  GOBP_WFC        = "GOBP\nWFC",
  KEGG_PPAR       = "KEGG\nPPAR"
)

draw_heatmap(heat_df, adipo_sig, adipo_annot, adipo_order,
             "Adipogenesis genes ~ BMI\ncommitted preadipocyte (nominal p<0.05)",
             file.path(OUT_DIR, "03.committed_aspc_adipogenesis_heatmap.pdf"),
             gs_labels = adipo_gs_labels,
             sig_stats = sig_stats %>% filter(gene %in% adipo_sig),
             tile_h    = 0.15)

cat("\nDone.\n")
