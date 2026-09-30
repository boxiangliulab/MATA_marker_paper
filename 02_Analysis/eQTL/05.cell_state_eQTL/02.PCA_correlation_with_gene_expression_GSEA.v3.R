suppressPackageStartupMessages({
  library(getopt)
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(fgsea)
  library(msigdbr)
  library(readr)
  library(purrr)
  library(ggplot2)
  library(ggrastr)
  library(broom)
  library(data.table)

  # for gseGO
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

## ============ CLI ============
command = matrix(c(
  "depot","d",1,"character",
  "celltype","c",1,"character",
  "help","h",0,"logical"
), byrow = TRUE, ncol = 4)
args = getopt(command)
if (!is.null(args$help) || is.null(args$depot)) {
  cat("Usage: Rscript script.R -d <Depot> [-c <Celltype>]\n"); quit(save="no", status=1)
}
depot <- args$depot      # "SAT" / "VAT"
# celltype <- args$celltype  # not used in this script version

# Project root is 8 levels up from this script
args0 <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args0[grep("--file=", args0)])
project_root <- normalizePath(file.path(dirname(normalizePath(script_path)),
                                        paste(rep("..", 8), collapse = "/")))

barcode_col <- "barcode"

## ============ Palette (your 11 colors) ============
gsea_colors <- c(
  "#a50026",
  "#d73027",
  "#f46d43",
  "#fdae61",
  "#fee090",
  "#ffffbf",
  "#e0f3f8",
  "#abd9e9",
  "#74add1",
  "#4575b4",
  "#313695"
)
# ensure blue=negative, red=positive for gradient scale
gsea_colors_blue2red <- rev(gsea_colors)  # "#313695" ... "#a50026"

## ============ Config ============
set.seed(123)
n_pc           <- 8
gs_source      <- "GO"
minSize        <- 10
maxSize        <- 500
topN_pathway   <- 20
topK_genes     <- 50
corr_method    <- "pearson"
min_frac_expr  <- 0.01
assay_for_cor  <- "SCT"

base_dir <- file.path(project_root, "data/09.eQTL/03.Saige_QTL/02.output/11.cell_state_eQTL/02.adipogenesis")
obj_path <- file.path(base_dir, depot)

out_root <- file.path(base_dir, depot, "output")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

gsea_outdir  <- file.path(out_root, "GSEA")
plot_outdir  <- file.path(out_root, "GSEA", "plot")
dir.create(gsea_outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_outdir, recursive = TRUE, showWarnings = FALSE)

# -----------------------
# Read 10x-like matrix (genes x cells)
# -----------------------
options(future.globals.maxSize = 500 * 1024^2)
mat <- ReadMtx(
  mtx      = file.path(obj_path, "counts.mtx"),
  features = file.path(obj_path, "genes.tsv"),
  cells    = file.path(obj_path, "barcodes.tsv")
)

obj <- CreateSeuratObject(counts = mat, assay = "RNA", project = paste0(depot, "_ASPC_and_adipocyte"))
rm(mat); gc()

message("Seurat object created: Cells=", ncol(obj), " Genes=", nrow(obj))

# -----------------------
# Merge metadata (per cell)
# -----------------------
metadata_tsv <- file.path(project_root, "data/08.integration_snRNA_object/01.freeze_version/metadata_after_libraryQC_after_ctQC_freeze.tsv")
meta <- fread(metadata_tsv)

if (!(barcode_col %in% colnames(meta))) {
  stop("Barcode column '", barcode_col, "' not found in metadata TSV. Columns include:\n",
       paste(head(colnames(meta), 50), collapse = ", "), "\n... (showing first 50)")
}

meta[[barcode_col]] <- as.character(meta[[barcode_col]])
meta2 <- meta %>% distinct(.data[[barcode_col]], .keep_all = TRUE)
meta2 <- meta2[meta2[[barcode_col]] %in% colnames(obj), ]
message("Metadata rows matched to object: ", nrow(meta2))
meta2 <- as.data.frame(meta2)
rownames(meta2) <- meta2[[barcode_col]]
obj <- AddMetaData(obj, metadata = as.data.frame(meta2))
obj <- subset(obj, Annotation_level1=="ASPC")
# -----------------------
# Add QC covariates used in SCT regression
# -----------------------
obj[["percent.mt"]]   <- PercentageFeatureSet(obj, pattern = "^MT-")
obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = "^[Rr][Pp][LlSs]")

# ============================================================
# Normalize + Cell cycle (keep this order)
# ============================================================
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 1e4, verbose = FALSE)

cc.genes <- Seurat::cc.genes.updated.2019
obj <- CellCycleScoring(
  obj,
  s.features = cc.genes$s.genes,
  g2m.features = cc.genes$g2m.genes,
  set.ident = FALSE
)

# ============================================================
# SCTransform + PCA
# ============================================================
options(future.globals.maxSize = 20000 * 1024^2)
obj <- SCTransform(
  obj,
  assay = "RNA",
  vars.to.regress = c("percent.mt", "percent.ribo", "S.Score", "G2M.Score"),
  verbose = FALSE
)
obj <- RunPCA(obj, assay = "SCT", npcs = n_pc, verbose = FALSE)

# Save PCA embeddings
pca_embed <- Embeddings(obj, reduction = "pca")[, 1:n_pc, drop = FALSE]
pca_out <- data.frame(barcode = rownames(pca_embed), pca_embed, check.names = FALSE)
pca_embed_fp <- file.path(out_root, sprintf("%s_%s_PCA_embedding_PC1-%d.tsv", depot, "adipogenesis", n_pc))
fwrite(pca_out, pca_embed_fp, sep = "\t")

# Save PCA loadings
loadings <- Loadings(obj, reduction = "pca")[, 1:n_pc, drop = FALSE]
fwrite(as.data.table(loadings, keep.rownames = "gene"),
       file.path(out_root, sprintf("%s_%s_PCA_loadings_PC1-%d.tsv", depot, "adipogenesis", n_pc)),
       sep = "\t")

# ============================================================
# Gene filtering for correlation (>= min_frac_expr nuclei expressed)
# ============================================================
rna_data <- GetAssayData(obj, assay = "RNA", layer = "data")
frac_expr <- Matrix::rowMeans(rna_data > 0)
keep_genes <- names(frac_expr)[frac_expr >= min_frac_expr]

message("Gene filter min_frac_expr=", min_frac_expr,
        " kept=", length(keep_genes), "/", length(frac_expr))

sct_scaled <- GetAssayData(obj, assay = "SCT", layer = "scale.data")
keep_genes <- intersect(keep_genes, rownames(sct_scaled))
sct_scaled <- sct_scaled[keep_genes, , drop = FALSE]

# ============================================================
# gseGO per PC
# ============================================================

run_one_pc_gseGO <- function(pc_name,
                             pc_vec,
                             expr_mat,
                             out_dir,
                             depot,
                             go_ont   = "BP",
                             go_minGS = 10,
                             go_maxGS = 800,
                             go_pcut  = 1) {

  common_cells <- intersect(names(pc_vec), colnames(expr_mat))
  pc_vec <- pc_vec[common_cells]
  expr_mat <- expr_mat[, common_cells, drop = FALSE]

  cor_res <- lapply(rownames(expr_mat), function(g) {
    x <- as.numeric(expr_mat[g, ])
    ok <- is.finite(x) & is.finite(pc_vec)
    if (sum(ok) < 10) return(c(cor = NA_real_, p = NA_real_))
    ct <- suppressWarnings(cor.test(x[ok], pc_vec[ok], method = "pearson"))
    c(cor = unname(ct$estimate), p = ct$p.value)
  })

  cor_dt <- data.table::rbindlist(lapply(seq_along(cor_res), function(i) {
    data.table::data.table(
      gene = rownames(expr_mat)[i],
      cor  = cor_res[[i]][["cor"]],
      p    = cor_res[[i]][["p"]]
    )
  }))

  data.table::fwrite(
    cor_dt,
    file.path(out_dir, sprintf("%s_%s_%s_cor_p.tsv", depot, "adipogenesis", pc_name)),
    sep = "\t"
  )

  #gene_list_go <- with(cor_dt, sign(cor) * (-log10(p)))
  gene_list_go <- cor_dt$cor
  names(gene_list_go) <- cor_dt$gene
  gene_list_go <- gene_list_go[!is.na(gene_list_go) & is.finite(gene_list_go)]
  gene_list_go <- gene_list_go[order(abs(gene_list_go), decreasing = TRUE)]
  gene_list_go <- gene_list_go[!duplicated(names(gene_list_go))]
  gene_list_go <- sort(gene_list_go, decreasing = TRUE)

  ego <- tryCatch(
    clusterProfiler::gseGO(
      geneList      = gene_list_go,
      ont           = go_ont,
      keyType       = "SYMBOL",
      OrgDb         = org.Hs.eg.db,
      minGSSize     = go_minGS,
      maxGSSize     = go_maxGS,
      pAdjustMethod = "fdr",
      pvalueCutoff  = go_pcut,
      verbose       = FALSE,
      seed          = TRUE,
      eps           = 0
    ),
    error = function(e) NULL
  )

  if (!is.null(ego) && nrow(as.data.frame(ego@result)) > 0) {
    res_go <- as.data.frame(ego@result)

    readr::write_csv(
      res_go,
      file.path(out_dir, sprintf("%s.%s.%s.gseGO_%s.csv", "adipogenesis", depot, pc_name, go_ont))
    )

    # ============================
    # NEW: Dotplot with NES (x=NES)
    # ============================
    dot_df <- res_go %>%
      dplyr::arrange(p.adjust) %>%
      dplyr::slice_head(n = topN_pathway) %>%
      dplyr::mutate(
        term = Description,
        logPadj = -log10(p.adjust)
      )

    max_abs <- max(abs(dot_df$NES), na.rm = TRUE)

    p_dot <- ggplot(dot_df, aes(x = NES, y = reorder(term, NES))) +
      geom_point(aes(size = logPadj, color = NES)) +
      scale_color_gradientn(
        colours = gsea_colors_blue2red,
        limits = c(-max_abs, max_abs),
        name = "NES"
      ) +
      scale_size_continuous(name = "-log10(adj.P)") +
      theme_classic(base_size = 12) +
      theme(axis.text.y = element_text(size = 9)) +
      labs(
        title = paste0("adipogenesis ", depot, " ", pc_name, " gseGO-", go_ont),
        x = "Normalized Enrichment Score (NES)",
        y = NULL
      )

    ggplot2::ggsave(
      file.path(out_dir, sprintf("%s.%s.%s.gseGO_%s.dotplot_NES.pdf",
                                 "adipogenesis", depot, pc_name, go_ont)),
      p_dot, width = 7.5, height = 6.5, dpi = 200
    )

    return(res_go)
  } else {
    message("[!] No significant gseGO terms for ", pc_name, " (", go_ont, ").")
    return(NULL)
  }
}

# ============================================================
# Run across PCs and save combined summary
# ============================================================
pc_mat <- Embeddings(obj, "pca")[, 1:n_pc, drop = FALSE]
all_summary <- list()

go_ont   <- "BP"
go_minGS <- 10
go_maxGS <- 800
go_pcut  <- 1

######  GSEA GO
for (k in seq_len(n_pc)) {
  pc_name <- paste0("PC", k)
  message("Running: ", pc_name)

  pc_vec <- pc_mat[, k]
  names(pc_vec) <- rownames(pc_mat)

  res_go <- run_one_pc_gseGO(
    pc_name  = pc_name,
    pc_vec   = pc_vec,
    expr_mat = sct_scaled,
    out_dir  = gsea_outdir,
    depot    = depot,
    go_ont   = go_ont,
    go_minGS = go_minGS,
    go_maxGS = go_maxGS,
    go_pcut  = go_pcut
  )

  if (!is.null(res_go) && nrow(res_go) > 0) {
    summ_go <- res_go %>%
      dplyr::transmute(
        PC         = pc_name,
        collection = paste0("gseGO_", go_ont),
        ID         = ID,
        term       = Description,
        NES        = NES,
        p.adjust   = p.adjust,
        setSize    = setSize
      )
    all_summary[[pc_name]] <- summ_go
  }
}

all_summary_dt <- dplyr::bind_rows(all_summary)

data.table::fwrite(
  all_summary_dt,
  file.path(out_root, sprintf("%s_%s_allPC_gseGO_%s_summary.tsv", depot, "adipogenesis", go_ont)),
  sep = "\t"
)

# ============================================================
# NEW: Tile overview (padj<0.05 only) + global sorting of terms
# ============================================================
if (nrow(all_summary_dt) > 0) {

  top_map_sig <- all_summary_dt %>%
    dplyr::filter(p.adjust < 0.05)

  if (nrow(top_map_sig) == 0) {
    message("[!] No significant terms (padj<0.05) across PCs; skip tile plot.")
  } else {

    # keep top N per PC, but only among significant
    top_map_sig <- top_map_sig %>%
      dplyr::group_by(PC, collection) %>%
      dplyr::arrange(p.adjust) %>%
      dplyr::slice_head(n = 10) %>%
      dplyr::ungroup()

    # global ordering per collection:
    # 1) best p.adjust across PCs (min)
    # 2) tie-break: larger |mean NES|
    term_order_tbl <- top_map_sig %>%
      dplyr::group_by(collection, term) %>%
      dplyr::summarise(
        best_padj = min(p.adjust, na.rm = TRUE),
        mean_NES  = mean(NES, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::arrange(collection, best_padj, dplyr::desc(abs(mean_NES)))

    top_map_sig <- top_map_sig %>%
      dplyr::left_join(term_order_tbl, by = c("collection", "term")) %>%
      dplyr::group_by(collection) %>%
      dplyr::mutate(
        term = factor(term, levels = unique(term_order_tbl$term[term_order_tbl$collection == unique(collection)]))
      ) %>%
      dplyr::ungroup()

    max_abs_all <- max(abs(top_map_sig$NES), na.rm = TRUE)

    p_tile <- ggplot(top_map_sig, aes(x = PC, y = term, fill = NES)) +
      geom_tile(color = "white", linewidth = 0.2) +
      facet_wrap(~collection, scales = "free_y", ncol = 1) +
      scale_fill_gradientn(
        colours = gsea_colors_blue2red,
        limits  = c(-max_abs_all, max_abs_all),
        name    = "NES"
      ) +
      theme_classic(base_size = 12) +
      theme(
        axis.text.y = element_text(size = 6),
        panel.grid  = element_blank()
      ) +
      labs(
        title = paste0("Significant GO terms per PC (padj<0.05; tile=NES; ont=", go_ont, ")"),
        x = NULL, y = NULL
      )

    ggsave(
      file.path(gsea_outdir, sprintf("%s_%s_gseGO_%s_tile_overview_padj0.05.sorted.pdf",
                                     depot, "adipogenesis", go_ont)),
      p_tile, width = 9, height = 10.5, dpi = 200
    )
  }

} else {
  message("[!] No gseGO results across PCs; skip tile plot.")
}

#===============================
# UMAP visualization (fix PCA file reading)
#===============================
meta_path <- metadata_tsv
if (depot == "SAT"){
  umap_path <- file.path(project_root, "data/09.eQTL/03.Saige_QTL/02.output/11.cell_state_eQTL/02.adipogenesis/SAT/ASPC_scvi_reumap_from_allobj/SAT_ASPC_scVI_reUMAP.coords.tsv")
} else if (depot=="VAT"){
  umap_path <- file.path(project_root, "data/09.eQTL/03.Saige_QTL/02.output/11.cell_state_eQTL/02.adipogenesis/VAT/ASPC_scvi_reumap_from_allobj/VAT_ASPC_scVI_reUMAP.coords.tsv")
}

if (file.exists(meta_path) && file.exists(umap_path) && file.exists(pca_embed_fp)) {
  meta <- read.table(meta_path, sep = "\t", header = TRUE, check.names = FALSE)
  scvi_umap <- read.table(umap_path, header = TRUE,sep="\t")
  names(scvi_umap) <- c("barcode","scvi_UMAP_1","scvi_UMAP_2")

  pca_df <- readr::read_tsv(pca_embed_fp, show_col_types = FALSE)
  # pca_df: barcode + PC_1..PC_15
  merge_pca_umap <- merge(pca_df, scvi_umap, by = "barcode")
  merge_pca_umap <- merge(merge_pca_umap, meta, by = "barcode")

  readr::write_tsv(merge_pca_umap,
                   file.path(plot_outdir, paste0(depot,"-","adipogenesis","-merge_pca_umap.tsv")))

  pdf(file.path(plot_outdir,"UMAP_by_annotation_level2.pdf"), width = 6, height = 4)
  p1 <- ggplot(merge_pca_umap, aes(scvi_UMAP_1, scvi_UMAP_2, color = Annotation_level2)) +
    ggrastr::geom_point_rast(size = 0.01) + theme_classic() +
    ggtitle(paste0(depot,": Annotation_level2 (", nrow(merge_pca_umap) ,") nuclei")) + scale_color_manual(values= c( "Early_preadipocyte"= "#c7e9c0","Committed_preadipocyte" = "#a1d99b"))
  print(p1)
  dev.off()

  pc_cols <- grep("^PC_\\d+$", names(merge_pca_umap), value = TRUE)
  pc_cols <- pc_cols[order(as.integer(sub("PC_", "", pc_cols)))]
  merge_pca_umap[pc_cols] <- lapply(merge_pca_umap[pc_cols], function(x) as.numeric(as.character(x)))

  df_long <- merge_pca_umap %>%
    dplyr::select(scvi_UMAP_1, scvi_UMAP_2, all_of(pc_cols)) %>%
    pivot_longer(all_of(pc_cols), names_to = "PC", values_to = "score") %>%
    filter(is.finite(score))

  rng <- max(abs(df_long$score), na.rm = TRUE)

  p_raw <- ggplot(df_long, aes(scvi_UMAP_1, scvi_UMAP_2, color = score)) +
    ggrastr::geom_point_rast(size = 0.05, alpha = 0.8) +
    coord_equal() +
    scale_colour_gradientn(
      colours = c("#313695","#ffffbf","#a50026"),
      limits = c(-rng, rng),
      name = "PC value"
    ) +
    facet_wrap(~ PC, ncol = 5) +
    labs(x = "scVI UMAP 1", y = "scVI UMAP 2",
         title = paste0(depot, "-", "adipogenesis", " (raw)")) +
    theme_classic(base_size = 10) +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          strip.text = element_text(face = "bold"))

  ggsave(file.path(plot_outdir, "UMAP_by_PC_raw.pdf"), p_raw, width = 14, height = 10)
    # ============================
# NEW: scaled PCA (-1 to 1) with custom palette
# ============================

df_scaled <- df_long %>%
  group_by(PC) %>%
  mutate(
    scale_factor = max(abs(score), na.rm = TRUE),
    score_scaled = ifelse(
      scale_factor > 0,
      pmax(pmin(score / scale_factor, 1), -1),
      0
    )
  ) %>%
  ungroup()

p_scaled <- ggplot(df_scaled, aes(scvi_UMAP_1, scvi_UMAP_2, color = score_scaled)) +
  ggrastr::geom_point_rast(size = 0.2, alpha = 0.8) +
  coord_equal() +
  scale_color_gradient2(
    limits   = c(-1, 1),
    midpoint = 0,
    low      = "blue",
    mid      = "white",
    high     = "red",
    name     = "PC (scaled -1..1)"
  ) +
  facet_wrap(~PC, ncol = 5) +
  labs(
    x = "scVI UMAP 1",
    y = "scVI UMAP 2",
    title = paste0(depot, "-", "adipogenesis", " (scaled -1..1)")
  ) +
  theme_classic(base_size = 10) +
  theme(
    axis.text  = element_blank(),
    axis.ticks = element_blank(),
    strip.text = element_text(face = "bold")
  )

ggsave(
  file.path(plot_outdir, "UMAP_by_PC_scaled_-1_1_customPalette.pdf"),
  p_scaled,
  width = 14,
  height = 10
)
}

message("Done. All GSEA outputs are under: ", gsea_outdir)

## =============== Save phenotype for cell-state QTL analysis ============
DefaultAssay(obj) <- "RNA"
#celltypes <- c("Adipocyte","Committed_preadipocyte","Early_preadipocyte")
celltypes <- c("Committed_preadipocyte","Early_preadipocyte")
sig_qtl <- read.table(
  file.path(project_root, "data/09.eQTL/03.Saige_QTL/02.output/03.formal_calling/ancestry_specific_version/All_eQTL_SAT_and_VAT_with_beta.tsv"),
  sep = "\t", header = TRUE
)
sig_qtl <- sig_qtl[which(sig_qtl$celltype %in% celltypes & sig_qtl$depot == depot),]
gene_list <- unique(sig_qtl$gene)

qtl_expre_mat <- Seurat::GetAssayData(obj, assay = "RNA", layer = "counts")
qtl_expre_mat <- t(qtl_expre_mat) %>% as.data.frame()
qtl_expre_mat <- qtl_expre_mat[, gene_list, drop = FALSE]
qtl_expre_mat <- rownames_to_column(qtl_expre_mat, "barcode")
gene_col <- setdiff(names(qtl_expre_mat), "barcode")

# Merge with PC (use existing PCA embeddings in object)
pca_mat <- Embeddings(obj, "pca")
if (ncol(pca_mat) < n_pc) stop("Requested PCs exceed available PCs.")
pca_mat_sub <- pca_mat[, 1:n_pc, drop = FALSE]

qtl_pca <- as.data.frame(pca_mat_sub)
qtl_pca <- rownames_to_column(qtl_pca, "barcode")
names(qtl_pca) <- gsub("^PC_(\\d+)$", "pPC\\1", names(qtl_pca))
pca_col <- setdiff(names(qtl_pca), "barcode")

# FIX: this line was broken in your paste
qtl_pheno <- merge(qtl_expre_mat, qtl_pca, by = "barcode")

# Merge donor-level covariates from metadata (already in obj@meta.data)
meta_s <- as.data.frame(obj@meta.data[, c("nCount_RNA","percent.mt","nFeature_RNA","Annotation_level2","donor_id","Age","Gender","Race","BMI","Self_defined_DM")])
meta_s$barcode <- rownames(meta_s)
meta_col <- names(meta_s)

qtl_pheno <- merge(meta_s, qtl_pheno, by = "barcode")

# Genotype PCs
genotype_pc <- read.csv("/home/project/11003054/e1353486/project/adipose/fullgenotype/Ancestry_specific/01.PCA/296_sata_pca_ethnicity_result.csv")

if (depot == "SAT") {
  genotype_pc <- genotype_pc[, c("FID", paste0("PC", 1:8))]
  colnames(genotype_pc) <- gsub("^PC(\\d+)", "gPC\\1", colnames(genotype_pc))
  gpc_col <- paste0("gPC", 1:8)
} else {
  genotype_pc <- genotype_pc[, c("FID", paste0("PC", 1:6))]
  colnames(genotype_pc) <- gsub("^PC(\\d+)", "gPC\\1", colnames(genotype_pc))
  gpc_col <- paste0("gPC", 1:6)
}

merge_dt <- merge(qtl_pheno, genotype_pc, by.x = "donor_id", by.y = "FID")

merge_dt <- merge_dt[, c(meta_col, pca_col, gpc_col, gene_col)]
names(merge_dt) <- gsub("-", "_", names(merge_dt))

write.table(
  merge_dt,
  file.path(out_root, "celltype_state_phenotype.RNA.safe_name.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

message("Phenotype written: ", file.path(out_root, "celltype_state_phenotype.RNA.safe_name.tsv"))
