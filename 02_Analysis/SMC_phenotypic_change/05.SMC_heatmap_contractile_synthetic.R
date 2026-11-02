suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(edgeR)
})

# ─── Config ───────────────────────────────────────────────────────────────────
# Project root is 7 levels up from this script
args0 <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args0[grep("--file=", args0)])
project_root <- normalizePath(file.path(dirname(normalizePath(script_path)),
                                        paste(rep("..", 7), collapse = "/")))

OUT_DIR <- file.path(project_root, "data/35.vascular/03.SMC_heatmap")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

METADATA_PATH  <- file.path(project_root, "data/08.integration_snRNA_object/01.freeze_version/metadata_after_libraryQC_after_ctQC_freeze.tsv")
BASE_COUNT_DIR <- file.path(project_root, "data/08.integration_snRNA_object/scvi/03.count_matrix")
COUNT_SAT <- file.path(BASE_COUNT_DIR, "SAT/01.raw_count/SMC.raw_count.tsv")
COUNT_VAT <- file.path(BASE_COUNT_DIR, "VAT/01.raw_count/SMC.raw_count.tsv")
GENE_DIR  <- file.path(project_root, "data/35.vascular/02.module_score")  # kept for reference

ASIAN_RACES  <- c("Chinese", "Malay", "Indian")
MIN_CELLS    <- 10
COL_VAT      <- "#F3AF44"
COL_SAT      <- "#67B3AD"
DEPOT_COLORS <- c(SAT = COL_SAT, VAT = COL_VAT)
GROUP_COLORS <- c(
  "Non-obese"      = "#efedf5",
  "Obese"          = "#bcbddc",
  "Severely obese" = "#807dba",
  "Morbidly obese" = "#54278f"
)
BMI_LEVELS <- c("Non-obese", "Obese", "Severely obese", "Morbidly obese")

# ─── Read gene lists ──────────────────────────────────────────────────────────
MODULE_SCORE_DIR <- file.path(project_root, "data/35.vascular/04.SMC_module_score")
SMOOTH_GENES <- fread(file.path(MODULE_SCORE_DIR,
                                "00a.GOBP_SMOOTH_MUSCLE_CONTRACTION_leading_genes.tsv"))$gene
ECM_GENES    <- fread(file.path(MODULE_SCORE_DIR,
                                "00b.GOBP_COLLAGEN_FIBRIL_ORGANIZATION_leading_genes.tsv"))$gene

# No overlap between the two sets
GENE_ORDER <- c(SMOOTH_GENES, setdiff(ECM_GENES, SMOOTH_GENES))
GENE_ANNOT <- setNames(
  ifelse(GENE_ORDER %in% SMOOTH_GENES & GENE_ORDER %in% ECM_GENES, "Both",
         ifelse(GENE_ORDER %in% SMOOTH_GENES, "Contractile", "Collagen/ECM")),
  GENE_ORDER
)
message("Display order (", length(GENE_ORDER), " genes): ", paste(GENE_ORDER, collapse = ", "))

ANNOT_COLORS <- c(Contractile = "#2171b5", Both = "#74c476", "Collagen/ECM" = "#e6550d")
HEAT_COLORS  <- c("#313695","#4575b4","#74add1","#abd9e9","#e0f3f8",
                  "#ffffbf","#fee090","#fdae61","#f46d43","#d73027","#a50026")

# ─── Load metadata ────────────────────────────────────────────────────────────
message("[1] Loading metadata...")
meta <- fread(METADATA_PATH, sep = "\t", header = TRUE, data.table = FALSE)
cell_meta <- meta %>%
  select(barcode, donor_id, Depot, Race, BMI, Age, Gender, Diabetic) %>%
  filter(Race %in% ASIAN_RACES, !is.na(BMI))

donor_meta <- cell_meta %>%
  distinct(donor_id, .keep_all = TRUE) %>%
  mutate(
    T2DM          = ifelse(Diabetic == "T2DM", 1, 0),
    BMI_group     = case_when(
      BMI <  27.5 ~ "Non-obese",
      BMI <  32.5 ~ "Obese",
      BMI <  37.5 ~ "Severely obese",
      TRUE        ~ "Morbidly obese"
    ),
    BMI_group     = factor(BMI_group, levels = BMI_LEVELS),
    BMI_group_num = as.numeric(BMI_group)
  )

# ─── Full-matrix pseudobulk (ALL genes) for proper edgeR TMM ──────────────────
# Standard pseudobulk: sum ALL gene counts per donor × depot,
# then edgeR filterByExpr + calcNormFactors(TMM) on the full matrix.
# This ensures TMM normalization is computed from the full transcriptome,
# not biased by a small target gene subset.
message("[2] Loading full count matrices and computing pseudobulk...")

pseudobulk_full <- function(fp, depot, cell_meta, donors_keep) {
  message("  Reading: ", fp)
  dt <- fread(fp, sep = "\t", data.table = TRUE)
  bc_col <- names(dt)[1]

  # Filter to valid Asian cells of retained donors
  valid_bc <- cell_meta %>%
    filter(Depot == depot, donor_id %in% donors_keep) %>%
    pull(barcode)
  dt <- dt[dt[[bc_col]] %in% valid_bc]

  # Add donor_id
  bc2donor <- setNames(
    cell_meta$donor_id[match(valid_bc, cell_meta$barcode)],
    valid_bc
  )
  dt[, donor_id := bc2donor[dt[[bc_col]]]]
  dt[, (bc_col) := NULL]

  # Pseudobulk: sum per donor across ALL genes
  pb <- dt[, lapply(.SD, sum, na.rm = TRUE), by = donor_id]
  donors <- pb$donor_id
  pb[, donor_id := NULL]
  mat <- t(as.matrix(pb))   # genes × donors
  colnames(mat) <- donors
  message("  Done: ", ncol(mat), " donors, ", nrow(mat), " genes")
  mat
}

# First determine donors_keep (need cell counts from both depots)
# Use a quick read of barcode + donor mapping from metadata
message("  Determining valid donors (≥", MIN_CELLS, " cells in both depots)...")
cnt_sat <- cell_meta %>% filter(Depot == "SAT") %>% count(donor_id, name = "n_sat")
cnt_vat <- cell_meta %>% filter(Depot == "VAT") %>% count(donor_id, name = "n_vat")
donors_keep <- inner_join(cnt_sat, cnt_vat, by = "donor_id") %>%
  filter(n_sat >= MIN_CELLS, n_vat >= MIN_CELLS) %>%
  pull(donor_id)
message("  Donors retained: ", length(donors_keep))

pb_sat <- pseudobulk_full(COUNT_SAT, "SAT", cell_meta, donors_keep)
pb_vat <- pseudobulk_full(COUNT_VAT, "VAT", cell_meta, donors_keep)

# Combine into one matrix: genes × (all donors × depots)
common_genes <- intersect(rownames(pb_sat), rownames(pb_vat))
colnames(pb_sat) <- paste(colnames(pb_sat), "SAT", sep = "_")
colnames(pb_vat) <- paste(colnames(pb_vat), "VAT", sep = "_")
pb_all <- cbind(pb_sat[common_genes, ], pb_vat[common_genes, ])

pb_info <- data.frame(
  pb_id    = colnames(pb_all),
  donor_id = c(gsub("_SAT$", "", colnames(pb_sat)),
               gsub("_VAT$", "", colnames(pb_vat))),
  Depot    = c(rep("SAT", ncol(pb_sat)), rep("VAT", ncol(pb_vat)))
)

# ─── edgeR: standard pipeline on FULL gene set ────────────────────────────────
message("[3] edgeR: filterByExpr + calcNormFactors(TMM) on full gene set...")
dge <- DGEList(counts = pb_all)

# filterByExpr: keep genes expressed in enough samples
design_filter <- model.matrix(~ pb_info$Depot)
keep <- filterByExpr(dge, design = design_filter)
message("  Genes after filterByExpr: ", sum(keep), " / ", nrow(dge))
dge <- dge[keep, , keep.lib.sizes = FALSE]

# TMM normalization on the full (filtered) gene set
dge <- calcNormFactors(dge, method = "TMM")
message("  TMM norm factors range: [",
        round(min(dge$samples$norm.factors), 3), ", ",
        round(max(dge$samples$norm.factors), 3), "]")

# log2(CPM+1) for all genes → extract target genes
log_cpm_full <- cpm(dge, log = TRUE, prior.count = 1)

# ─── Extract target genes ─────────────────────────────────────────────────────
gene_cols <- intersect(GENE_ORDER, rownames(log_cpm_full))
missing_g <- setdiff(GENE_ORDER, rownames(log_cpm_full))
if (length(missing_g) > 0)
  message("  Genes not passing filterByExpr: ", paste(missing_g, collapse = ", "))

log_cpm <- log_cpm_full[gene_cols, , drop = FALSE]

# ─── Per-donor expression table ───────────────────────────────────────────────
expr_long <- as.data.frame(t(log_cpm)) %>%
  mutate(pb_id = rownames(.)) %>%
  left_join(pb_info, by = "pb_id") %>%
  left_join(donor_meta %>%
              select(donor_id, BMI, BMI_group, BMI_group_num, Age, Gender, Race, T2DM),
            by = "donor_id") %>%
  filter(!is.na(BMI_group)) %>%
  mutate(Depot     = factor(Depot, levels = c("SAT","VAT")),
         BMI_group = factor(BMI_group, levels = BMI_LEVELS))

write.table(expr_long,
            file.path(OUT_DIR, "01b.SMC_contractile_synthetic_logCPM_per_donor.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ─── LM per gene × depot ──────────────────────────────────────────────────────
# Model: log2(TMM-CPM+1) ~ BMI_group_num + Age + Gender + T2DM + Race
message("[4] LM significance tests...")
safe_lm <- function(df) {
  covars <- c("Age","Gender","T2DM","Race")
  keep_cov <- covars[sapply(covars, function(v) {
    x <- df[[v]]
    if (is.numeric(x)) sd(x, na.rm=TRUE) > 0 else length(unique(na.omit(x))) >= 2
  })]
  fmla <- as.formula(paste("log_cpm ~ BMI_group_num",
    if (length(keep_cov) > 0) paste("+", paste(keep_cov, collapse=" + ")) else ""))
  tryCatch({
    m <- lm(fmla, data=df)
    ct <- summary(m)$coefficients
    data.frame(n=nrow(df),
               beta  = ct["BMI_group_num","Estimate"],
               se    = ct["BMI_group_num","Std. Error"],
               t_val = ct["BMI_group_num","t value"],
               p_val = ct["BMI_group_num","Pr(>|t|)"])
  }, error=function(e) data.frame(n=nrow(df), beta=NA, se=NA, t_val=NA, p_val=NA))
}

sig_df <- expr_long %>%
  pivot_longer(all_of(gene_cols), names_to="gene", values_to="log_cpm") %>%
  group_by(gene, Depot) %>%
  group_modify(~safe_lm(.x)) %>%
  ungroup() %>%
  mutate(
    p_adj  = p.adjust(p_val, method="BH"),
    stars  = case_when(
      is.na(p_adj)  ~ "?",
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE          ~ "\u00d7"
    ),
    gene  = factor(gene, levels=rev(gene_cols)),
    Depot = factor(Depot, levels=c("SAT","VAT"))
  )

write.table(sig_df,
            file.path(OUT_DIR, "01d.SMC_lm_significance.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
message("  Results (p<0.1):")
print(sig_df %>% filter(p_val < 0.1) %>%
        select(gene, Depot, beta, p_val, p_adj, stars) %>%
        arrange(p_val) %>% as.data.frame())

# ─── Z-score for visualization only ──────────────────────────────────────────
expr_pivot <- expr_long %>%
  pivot_longer(all_of(gene_cols), names_to="gene", values_to="log_cpm")

expr_z <- expr_pivot %>%
  group_by(gene, Depot) %>%
  mutate(z_score = as.numeric(scale(log_cpm))) %>%
  ungroup()

heat_df <- expr_z %>%
  group_by(gene, Depot, BMI_group) %>%
  summarise(mean_z=mean(z_score, na.rm=TRUE), n_donors=n(), .groups="drop") %>%
  mutate(
    gene      = factor(gene, levels=rev(gene_cols)),
    BMI_group = factor(BMI_group, levels=BMI_LEVELS),
    Depot     = factor(Depot, levels=c("SAT","VAT"))
  )

write.table(heat_df,
            file.path(OUT_DIR, "01c.SMC_heatmap_zscored_data.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
n_per_group <- expr_long %>% count(Depot, BMI_group, name="n_donors")
write.table(n_per_group, file.path(OUT_DIR, "01e.SMC_n_donors_per_group.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

# ─── Plot: Heatmap ────────────────────────────────────────────────────────────
message("[5] Plotting heatmap...")
n_genes <- length(gene_cols)
n_bmi   <- 4
n_ecm   <- sum(gene_cols %in% setdiff(ECM_GENES, SMOOTH_GENES))
sep_y   <- n_ecm + 0.5   # ECM genes are at the bottom after rev(), separator above them

sig_plot_df <- sig_df %>% mutate(is_ns = p_adj >= 0.05 | is.na(p_adj))

# Gene-type annotation sidebar
annot_df <- data.frame(
  gene  = factor(gene_cols, levels=rev(gene_cols)),
  annot = factor(GENE_ANNOT[gene_cols], levels=c("Contractile","Both","Collagen/ECM"))
)
p_annot <- ggplot(annot_df, aes(x=0.5, y=gene, fill=annot)) +
  geom_tile(color="white", linewidth=0.4, width=1, height=1) +
  scale_fill_manual(values=ANNOT_COLORS, name=NULL) +
  scale_x_continuous(expand=c(0,0), limits=c(0,1)) +
  labs(x=NULL, y=NULL) +
  theme_void(base_size=8) +
  theme(
    legend.position  = "none",
    axis.text.y      = element_text(size=7, hjust=1, margin=margin(r=1)),
    plot.margin      = margin(0, 0, 0, 0)
  )

# Heatmap tiles
p_heat <- ggplot(heat_df, aes(x=BMI_group, y=gene, fill=mean_z)) +
  geom_tile(color="white", linewidth=0.35) +
  geom_hline(yintercept=sep_y, color="grey35", linewidth=0.6, linetype="dashed") +
  facet_wrap(~Depot, nrow=1) +
  scale_fill_gradientn(colors=HEAT_COLORS, name="Mean\nz-score", na.value="grey90") +
  labs(x=NULL, y=NULL,
       title="SMC contractile / collagen-ECM genes (GSEA leading edge)") +
  theme_classic(base_size=8) +
  theme(
    axis.text.x       = element_text(size=7, angle=30, vjust=1, hjust=1),
    axis.text.y       = element_blank(),
    axis.ticks.y      = element_blank(),
    axis.line         = element_blank(),
    strip.background  = element_blank(),
    strip.text        = element_text(face="bold", size=8.5),
    legend.key.height = unit(0.6,"cm"),
    legend.key.width  = unit(0.25,"cm"),
    legend.text       = element_text(size=7),
    legend.title      = element_text(size=7, face="bold"),
    plot.title        = element_text(size=8, face="bold"),
    plot.margin       = margin(0, 0, 0, 0)
  )

# Significance panel — square tiles
p_sig <- ggplot(sig_plot_df, aes(x=Depot, y=gene)) +
  geom_tile(fill="white", color="grey75", linewidth=0.3, width=1, height=1) +
  geom_text(data=sig_plot_df %>% filter( is_ns),
            aes(label=stars), color="grey65", size=2.5) +
  geom_text(data=sig_plot_df %>% filter(!is_ns, Depot=="SAT"),
            aes(label=stars), color=COL_SAT, size=2.8, fontface="bold") +
  geom_text(data=sig_plot_df %>% filter(!is_ns, Depot=="VAT"),
            aes(label=stars), color=COL_VAT, size=2.8, fontface="bold") +
  scale_x_discrete(expand=c(0,0)) +
  labs(x=NULL, y=NULL, title="LM sig.\n(BH-adj. p)") +
  theme_void(base_size=7) +
  theme(
    axis.text.x     = element_text(size=6.5, angle=30, vjust=1, hjust=1,
                                   color=c(COL_SAT, COL_VAT)),
    axis.text.y     = element_blank(),
    plot.title      = element_text(size=7, face="bold", hjust=0.5),
    plot.margin     = margin(0, 0, 0, 1)
  )

# Phenotype legend panel (separate, placed left of annotation)
annot_legend_df <- data.frame(
  x=1, y=c(3,2,1),
  label=c("Contractile","Both","Collagen/ECM"),
  fill=c("#2171b5","#74c476","#e6550d")
)
p_legend <- ggplot(annot_legend_df, aes(x=x, y=y, fill=label)) +
  geom_tile(width=0.6, height=0.6, color="white") +
  geom_text(aes(x=1.45, label=label), hjust=0, size=2.3, color="grey20") +
  scale_fill_manual(values=ANNOT_COLORS) +
  scale_x_continuous(limits=c(0.5, 3)) +
  coord_equal() +
  theme_void() +
  theme(legend.position="none", plot.margin=margin(0,0,0,0))

# Figure caption
cap_heat <- paste0(
  "Color = mean z-score of pseudobulk log\u2082(TMM-CPM+1) across donors ",
  "(z-scored per gene per depot). x = BMI group (1=Non-obese BMI<27.5, ",
  "2=Obese 27.5\u201332.5, 3=Severely obese 32.5\u201337.5, 4=Morbidly obese \u226537.5). ",
  "Right panel: nominal LM p-value (uncorrected) for BMI_group_num term ",
  "(log\u2082CPM ~ BMI_group_num + Age + Gender + T2DM + Race; ",
  "pseudobulk TMM-normalized on full transcriptome; n=", length(donors_keep), " donors). ",
  "Stars: * p<0.05, ** p<0.01, *** p<0.001; \u00d7 = not significant. ",
  "Star color = depot (SAT=teal, VAT=orange)."
)

# Combine: annotation | heatmap | sig panel
# Widths proportional to number of tiles: 1 : (4*2) : 2
p_combined <- (p_annot + p_heat + p_sig) +
  plot_layout(widths=c(1, n_bmi*2, 2), guides="keep") +
  plot_annotation(
    caption = cap_heat,
    theme   = theme(
      plot.caption = element_text(size=5.5, color="grey40", hjust=0, lineheight=1.3),
      plot.margin  = margin(2, 2, 2, 2)
    )
  ) &
  theme(plot.margin = margin(0, 0, 0, 0))

# Tile height ~0.10 in (half of previous), width kept proportional
tile_w <- 0.20   # column width
tile_h <- 0.10   # row height (narrow)
fig_h  <- n_genes * tile_h + 1.6
fig_w  <- (1 + n_bmi*2 + 2) * tile_w + 3.2

ggsave(file.path(OUT_DIR, "01.SMC_contractile_synthetic_heatmap.pdf"),
       p_combined, width=fig_w, height=fig_h)
message("  Saved heatmap: ", round(fig_w,1), " × ", round(fig_h,1), " in")

# ─── Plot: Violin plots for LM-significant genes ──────────────────────────────
message("[6] Violin plots for significant genes...")
sig_genes_any <- sig_df %>% filter(p_adj < 0.05) %>%
  pull(gene) %>% as.character() %>% unique()
message("  Significant genes: ", paste(sig_genes_any, collapse=", "))

trend_info <- sig_df %>%
  mutate(trend = ifelse(!is.na(p_adj) & p_adj < 0.05, "sig", "ns"),
         gene  = as.character(gene))

lm_labels <- sig_df %>%
  mutate(gene  = as.character(gene),
         label = paste0("adj.p=", formatC(p_adj, format="e", digits=1)))

n_by_group <- expr_long %>% count(Depot, BMI_group, name="n")

vln_plots <- lapply(sig_genes_any, function(g) {
  df_g  <- expr_pivot %>% filter(gene==g) %>%
    left_join(trend_info %>% filter(gene==g), by=c("gene","Depot"))
  lab_g <- lm_labels %>% filter(gene==g)

  ggplot(df_g, aes(x=BMI_group, y=log_cpm)) +
    geom_violin(aes(fill=Depot), trim=TRUE, scale="width",
                linewidth=0.25, alpha=0.5) +
    geom_jitter(aes(color=Depot), width=0.12, size=0.35, alpha=0.45) +
    geom_smooth(data=df_g %>% filter(trend=="sig"),
                aes(x=BMI_group_num, y=log_cpm),
                method="lm", formula=y~x, color="firebrick", fill="firebrick",
                linewidth=0.55, alpha=0.15, inherit.aes=FALSE) +
    geom_smooth(data=df_g %>% filter(trend=="ns"),
                aes(x=BMI_group_num, y=log_cpm),
                method="lm", formula=y~x, color="grey55", fill="grey70",
                linetype="dashed", linewidth=0.45, alpha=0.12, inherit.aes=FALSE) +
    geom_text(data=n_by_group,
              aes(x=BMI_group, y=-Inf, label=paste0("n=",n)),
              inherit.aes=FALSE, vjust=-0.3, size=1.8, color="grey45") +
    geom_text(data=lab_g,
              aes(x=2.5, y=Inf, label=label),
              inherit.aes=FALSE, hjust=0.5, vjust=1.2, size=1.9, lineheight=0.85) +
    facet_wrap(~Depot, nrow=1) +
    scale_fill_manual(values=DEPOT_COLORS, guide="none") +
    scale_color_manual(values=DEPOT_COLORS, guide="none") +
    labs(x=NULL, y=expression(log[2](TMM-CPM+1)), title=g) +
    theme_classic(base_size=8) +
    theme(axis.text.x     = element_text(size=6.5, angle=30, vjust=1, hjust=1),
          strip.background = element_blank(),
          strip.text       = element_text(face="bold", size=8),
          plot.title       = element_text(size=8.5, face="bold", hjust=0.5),
          plot.margin      = margin(3,3,3,3))
})
names(vln_plots) <- sig_genes_any

n_sig <- length(sig_genes_any)
ncols_vln <- min(n_sig, 4)

cap_vln <- paste0(
  "Each point = one donor's pseudobulk log\u2082(TMM-CPM+1) value (edgeR TMM-normalized ",
  "on full transcriptome). Violin = distribution across donors; jitter = individual donors. ",
  "Trend line: LM fit (red = p<0.05, grey dashed = p\u22650.05). ",
  "LM model: log\u2082CPM ~ BMI_group_num + Age + Gender + T2DM + Race. ",
  "Label: beta (effect per unit BMI group), stars and nominal p-value. n per group shown below x-axis."
)

p_vln <- wrap_plots(vln_plots, ncol=ncols_vln) +
  plot_annotation(
    title   = "SMC: LM-significant genes by BMI group",
    caption = cap_vln,
    theme   = theme(plot.title   = element_text(size=9, face="bold"),
                    plot.caption = element_text(size=5.5, color="grey40",
                                               hjust=0, lineheight=1.3))
  )

ggsave(file.path(OUT_DIR, "02.SMC_sig_genes_violin.pdf"),
       p_vln,
       width  = ncols_vln * 3.2,
       height = ceiling(n_sig / ncols_vln) * 3.0 + 1.0)
message("  Saved: 02.SMC_sig_genes_violin.pdf")
message("[Done] Script 01 complete.")
