
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(getopt)
  library(ggplot2)
  library(limma)

  # GSEA
  library(msigdbr)
  library(fgsea)
})

# =======================
# Args
# =======================
command <- matrix(c(
  "celltype",  "c", 1, "character",
  "help",      "h", 0, "logical"
), byrow = TRUE, ncol = 4)

args <- getopt(command)

if (!is.null(args$help) || is.null(args$celltype)) {
  cat(paste(getopt(command, usage = TRUE), "\n"))
  q(status = 1)
}

celltype <- args$celltype
# ============================================================
# Step 0) Settings / Paths
# ============================================================

metadata_path <- "../../../../../../../data/08.integration_snRNA_object/01.freeze_version/metadata_after_libraryQC_after_ctQC_freeze.tsv"

# Your method1 output (already exists)
method1_path <- paste0(
  "../../../../../../../data/21.differential_expressed_genes/01.pseudobulk/04.only_VAT_vs_SAT_pseudobulk/linear_model/flip_genes/",
  celltype, ".flip_genes.method1.tsv"
)

# CPM matrices (pseudobulk CPM) - you provided these
non_cpm_path <- paste0(
  "~/scratch/project/00.adipose/script/01.adipose_sqtl/21.differential_expressed_genes/by_comparisons/depot_difference_by_celltype_by_BMI_group/non_obese/non_obese.",
  celltype, ".tsv"
)
ob_cpm_path <- paste0(
  "~/scratch/project/00.adipose/script/01.adipose_sqtl/21.differential_expressed_genes/by_comparisons/depot_difference_by_celltype_by_BMI_group/obese/obese.",
  celltype, ".tsv"
)

# Output directory
out_dir <- "../../../../../../../data/21.differential_expressed_genes/01.pseudobulk/04.only_VAT_vs_SAT_pseudobulk/linear_model/interation_term/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# thresholds
bmi_cutoff <- 27.5
fdr_interaction_cutoff <- 0.05
topN_radar <- 20

# GSEA settings
gsea_nperm <- 20000
gsea_minSize <- 10
gsea_maxSize <- 500

# ============================================================
# Helpers
# ============================================================

read_gene_by_sample_tsv <- function(path) {
  dt <- fread(path)
  gene <- dt[[1]]
  m <- as.matrix(dt[, -1, with = FALSE])
  rownames(m) <- gene
  colnames(m) <- names(dt)[-1]
  storage.mode(m) <- "numeric"
  m
}
# parse sample_id to donor_id + Depot
# expected sample_id ends with _SAT or _VAT
parse_sample_id <- function(sample_ids) {
  depot <- sub("^.*_(SAT|VAT)$", "\\1", sample_ids)
  donor <- sub("_(SAT|VAT)$", "", sample_ids)
  data.frame(sample_id = sample_ids, donor_id = donor, Depot = depot, stringsAsFactors = FALSE)
}

safe_factor <- function(x) {
  x <- as.character(x)
  factor(x)
}

# ============================================================
# Step 1) Load donor metadata (from cell-level metadata, collapse to donor)
# ============================================================

message("[1/6] Loading metadata and building donor table ...")
metadata <- fread(metadata_path, sep = "\t", header = TRUE, data.table = FALSE)

metadata$BMI_group <- ifelse(metadata$BMI < bmi_cutoff, "non_obese", "obese")

meta_donor <- metadata %>%
  dplyr::select(donor_id, Age, BMI, Race, Gender, Diabetic, BMI_group) %>%
  dplyr::distinct(donor_id, .keep_all = TRUE)

# ============================================================
# Step 2) Load method1 table: flip_genes
# ============================================================

message("[2/6] Loading method1 gene table ...")
#if (!file.exists(method1_path)) stop("method1_path not found: ", method1_path)
#m1 <- fread(method1_path)

# expected columns (as you showed)
#req_cols <- c("gene","logFC_non","logFC_ob","delta_logFC")
#missing_cols <- setdiff(req_cols, colnames(m1))
#if (length(missing_cols) > 0) stop("method1 file missing columns: ", paste(missing_cols, collapse=", "))

#genes_method1 <- unique(m1$gene)

# ============================================================
# Step 3) Load CPM, merge, build sample meta, log2(CPM+1)
# ============================================================

message("[3/6] Loading CPM matrices ...")
#if (!file.exists(non_cpm_path)) stop("non_cpm_path not found: ", non_cpm_path)
#if (!file.exists(ob_cpm_path))  stop("ob_cpm_path not found: ",  ob_cpm_path)
if (file.exists(ob_cpm_path) & file.exists(non_cpm_path)) {
cpm_non <- read_gene_by_sample_tsv(non_cpm_path)
cpm_ob  <- read_gene_by_sample_tsv(ob_cpm_path)

# intersect genes, then cbind samples
common_genes <- intersect(rownames(cpm_non), rownames(cpm_ob))
if (length(common_genes) < 100) {
  warning("Common genes between non/obese CPM are small (", length(common_genes), "). Check inputs.")
}
cpm_non <- cpm_non[common_genes, , drop = FALSE]
cpm_ob  <- cpm_ob[common_genes,  , drop = FALSE]

cpm_all <- cbind(cpm_non, cpm_ob)
} else if (file.exists(ob_cpm_path) & !file.exists(non_cpm_path)){
        cpm_all <- read_gene_by_sample_tsv(ob_cpm_path)
        }


# sample-level meta parsed from colnames
sample_meta <- parse_sample_id(colnames(cpm_all))

# join donor covariates (BMI_group etc.)
sample_meta <- sample_meta %>%
  left_join(meta_donor, by = "donor_id")

if (any(is.na(sample_meta$BMI_group))) {
  bad <- sample_meta$donor_id[is.na(sample_meta$BMI_group)]
  stop("Some donors in CPM not found in donor meta. Example donor_id: ", bad[1])
}
# enforce factors
sample_meta$Depot     <- factor(sample_meta$Depot, levels = c("SAT","VAT"))
sample_meta$BMI_group <- factor(sample_meta$BMI_group, levels = c("non_obese","obese"))
sample_meta$donor_id  <- factor(sample_meta$donor_id)

# expression = log2(CPM+1)
expr <- log2(cpm_all + 1)

# filter genes tested
non_path <- file.path(
  "../../../../../../../data/21.differential_expressed_genes/01.pseudobulk/04.only_VAT_vs_SAT_pseudobulk/BMIgrp_non_obese/01.DEG/01.table",
  paste0(celltype, ".non_obese.paired.edgeR.all.tsv")
)

ob_path <- file.path(
  "../../../../../../../data/21.differential_expressed_genes/01.pseudobulk/04.only_VAT_vs_SAT_pseudobulk/BMIgrp_obese/01.DEG/01.table",
  paste0(celltype, ".obese.paired.edgeR.all.tsv")
)

exists_non <- file.exists(non_path)
exists_ob  <- file.exists(ob_path)

message("non_obese file: ", exists_non, " | ", non_path)
message("obese file:     ", exists_ob,  " | ", ob_path)

if (exists_non && exists_ob) {
  non_all <- read.table(non_path, sep = "\t", header = TRUE, check.names = FALSE)
  ob_all  <- read.table(ob_path,  sep = "\t", header = TRUE, check.names = FALSE)
  genes <- union(unique(ob_all$gene), unique(non_all$gene))
} else {
         ob_all  <- read.table(ob_path,  sep = "\t", header = TRUE, check.names = FALSE)
         genes <- unique(ob_all$gene)
}

expr <- expr[genes,]
# align meta to expr columns
stopifnot(all(sample_meta$sample_id == colnames(expr)))

# ============================================================
# Step 4) Linear model with interaction term (Depot * BMI), paired by donor
# ============================================================

message("[4/6] Fitting limma model for Depot * BMI interaction ...")

# Design: includes donor_id to handle pairing robustly (fixed effect)
design <- model.matrix(~ donor_id + Depot * BMI, data = sample_meta)
colnames(design) <- make.names(colnames(design))

# coef name for interaction
# After make.names, "DepotVAT:BMI" may become "DepotVAT.BMI"
coef_candidates <- grep("DepotVAT.*BMI", colnames(design), value = TRUE)
if (length(coef_candidates) != 1) {
  stop("Cannot uniquely find interaction coef. Found: ", paste(coef_candidates, collapse=", "))
}
coef_int <- coef_candidates[1]
message("Interaction coefficient = ", coef_int)

fit <- lmFit(expr, design)
fit <- eBayes(fit)

tt <- topTable(fit, coef = coef_int, number = Inf, sort.by = "none") %>%
  tibble::rownames_to_column("gene") %>%
  dplyr::mutate(
    Celltype = celltype,
    coef_tested = coef_int
  )

lm_out <- tt %>%
  dplyr::select(gene, Celltype, coef_tested,
                logFC, AveExpr, t, P.Value, adj.P.Val, B) %>%
  dplyr::rename(
    beta = logFC,
    p    = P.Value,
    fdr  = adj.P.Val
  )

# ============================================================
# Step 5) Intersect with method1 genes, pick top delta_logFC genes, radar plot
# ============================================================

message("[5/6] Selecting genes for radar plot ...")

lm_sig <- lm_out %>%
  filter(fdr < fdr_interaction_cutoff)
lm_out_file <- file.path(out_dir,"tables", paste0(celltype, ".DepotBMI.sig.limma.tsv"))
fwrite(lm_sig, lm_out_file, sep = "\t")

message("Saved interaction LM results: ", lm_out_file)
# intersect: (lm interaction significant) AND (in method1 list)
#genes_keep <- intersect(lm_sig$gene, genes_method1)

message("[5/6] Selecting genes for radar plot ...")


#if (length(genes_keep) == 0) {
#  stop("No genes found in intersection of (lm interaction significant) and method1 gene list. ",
#       "Try relaxing fdr_interaction_cutoff or check gene ID consistency.")
#}

#m1_keep <- m1 %>%
#  filter(gene %in% genes_keep) %>%
#  mutate(abs_delta = abs(delta_logFC)) %>%
#  arrange(desc(abs_delta)) %>%
#  slice_head(n = topN_radar)

#message("Radar genes selected: ", nrow(m1_keep))

# build radar df
#radar_df <- m1_keep %>%
#  dplyr::select(gene, logFC_non, logFC_ob) %>%
#  tidyr::pivot_longer(cols = c(logFC_non, logFC_ob),
#                      names_to = "group", values_to = "logFC") %>%
#  mutate(group = recode(group, logFC_non = "non_obese", logFC_ob = "obese"))

# order genes around circle by abs(delta) descending (or you can order by delta)
#radar_gene_levels <- m1_keep$gene
#radar_df$gene <- factor(radar_df$gene, levels = radar_gene_levels)

# To draw a clean radar with negative values, shift radius to be non-negative
#min_val <- min(radar_df$logFC, na.rm = TRUE)
#offset <- ifelse(min_val < 0, abs(min_val) + 0.5, 0)
#radar_df$r <- radar_df$logFC + offset

# reference circle at logFC=0 => r=offset
#ref_df <- data.frame(gene = factor(radar_gene_levels, levels = radar_gene_levels),
#                     r = rep(offset, length(radar_gene_levels)))

# y breaks (original scale)
#y_breaks <- pretty(radar_df$logFC, n = 5)
#y_breaks_r <- y_breaks + offset


#p_radar <- ggplot() +
#  geom_path(data = ref_df, aes(x = gene, y = r, group = 1), color = "red", linewidth = 0.6) +
#  geom_path(data = radar_df, aes(x = gene, y = r, group = group, color = group), linewidth = 0.9) +
#  geom_point(data = radar_df, aes(x = gene, y = r, color = group), size = 1.6) +
#  coord_polar() +
#  scale_y_continuous(breaks = y_breaks_r, labels = y_breaks) +
#  scale_color_manual(values = c(non_obese = "#8ecae6", obese = "#1d4ed8")) +
#  labs(
#    title = paste0(celltype, ": VAT vs SAT logFC in non-obese vs obese (top ", topN_radar, " by |delta_logFC|)"),
#    y = "logFC (VAT - SAT)",
#    x = NULL,
#    color = NULL
#  ) +
#  theme_minimal(base_size = 10) +
#  theme(
#    panel.grid.minor = element_blank(),
#    axis.text.x = element_text(size = 7),
#    plot.title = element_text(hjust = 0.5, face = "bold"),
#    legend.position = "top"
#  )

#radar_png <- file.path(out_dir, paste0(celltype, ".top", topN_radar, ".delta_logFC.radar.png"))
#ggsave(radar_png, p_radar, width = 8, height = 8, dpi = 200)
#message("Saved radar plot: ", radar_png)

# Also output the radar gene list used
#radar_gene_file <- file.path(out_dir, paste0(celltype, ".top", topN_radar, ".radar_genes.tsv"))
#fwrite(m1_keep, radar_gene_file, sep = "\t")
#message("Saved radar gene table: ", radar_gene_file)

# ============================================================
# Step 6) GSEA on interaction term (rank by t-stat from limma)
# ============================================================

message("[6/6] Running GSEA on interaction term ...")

# ranking by t-statistic is standard for preranked GSEA
tstat <- tt$t
names(tstat) <-tt$gene
tstat <- tstat[is.finite(tstat)]
tstat <- sort(tstat, decreasing = TRUE)
stat <- tstat

run_fgsea_save <- function(category, subcategory = NULL, tag, stat) {
  m <- msigdbr(species = "Homo sapiens", category = category)
  if (!is.null(subcategory)) m <- m %>% dplyr::filter(gs_subcat == subcategory)
  pathways <- split(m$gene_symbol, m$gs_name)

  fg <- fgsea::fgsea(
    pathways = pathways,
    stats    = stat,
    minSize  = gsea_minSize,
    maxSize  = gsea_maxSize,
    eps      = 1e-10   
  ) %>%
    as.data.frame() %>%
    dplyr::arrange(padj, dplyr::desc(abs(NES)))

  fg_sig <- fg %>% dplyr::filter(padj <= 0.05)

  out_tsv <- file.path(out_dir, "tables","GSEA", paste0(celltype, ".GSEA.sig.", tag, ".tsv"))
  data.table::fwrite(fg_sig, out_tsv, sep = "\t")

  list(res_all = fg, res_sig = fg_sig, out_tsv = out_tsv)
}


g_h <- run_fgsea_save("H", NULL, "Hallmark_interaction",stat)
g_bp <- run_fgsea_save("C5", "GO:BP", "GOBP_interaction",stat)

# Plot top pathways (5 positive + 5 negative by NES among padj<0.05)
plot_top_fgsea <- function(fg, tag) {
  fg_sig <- fg %>% filter(!is.na(padj)) %>% arrange(padj)
  fg_sig <- fg_sig %>% filter(padj < 0.05)

  if (nrow(fg_sig) == 0) {
    warning("No significant pathways for ", tag, " (padj<0.05). Skipping plots.")
    return(NULL)
  }

  pos <- fg_sig %>% filter(NES > 0) %>% arrange(padj) %>% slice_head(n = 5)
  neg <- fg_sig %>% filter(NES < 0) %>% arrange(padj) %>% slice_head(n = 5)
  pick <- bind_rows(pos, neg)

  pdf_file <- file.path(out_dir,"GSEA", paste0(celltype, ".GSEA.", tag, ".top_enrichment_plots.pdf"))
  pdf(pdf_file, width = 7, height = 5)

  for (pw in pick$pathway) {
    p <- plotEnrichment(msigdbr(species="Homo sapiens", category=ifelse(tag=="Hallmark","H","C5"),
                                subcategory=ifelse(tag=="Hallmark", NA, "GO:BP")) %>%
                          { if (tag=="Hallmark") . else filter(., gs_subcat=="GO:BP") } %>%
                          filter(gs_name == pw) %>%
                          pull(gene_symbol),
                        tstat) +
      ggtitle(paste0(pw, "\nNES=", round(pick$NES[pick$pathway==pw], 2),
                     ", padj=", signif(pick$padj[pick$pathway==pw], 3))) +
      theme_minimal(base_size = 11)
    print(p)
  }

  dev.off()
  message("Saved GSEA plots: ", pdf_file)
}

plot_top_fgsea(g_h$res, "Hallmark")
plot_top_fgsea(g_bp$res, "GOBP")

message("All done.")

