suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(Seurat)
  library(edgeR)
  library(getopt)
  library(pheatmap)
  library(reshape2)
  library(tibble)
  library(tidyverse)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(msigdbr)
  library(ggrepel)
  library(forcats)
  library(Matrix)
})

# =======================
# Args
# =======================
command <- matrix(c(
  "celltype",  "c", 1, "character",
  "bmi_group", "g", 1, "character",
  "help",      "h", 0, "logical"
), byrow = TRUE, ncol = 4)

args <- getopt(command)

if (!is.null(args$help) || is.null(args$celltype)  || is.null(args$bmi_group)) {
  cat(paste(getopt(command, usage = TRUE), "\n"))
  q(status = 1)
}

celltype  <- args$celltype
depot_tag <- args$depot
bmi_group <- args$bmi_group

# =======================
# Config
# =======================
min_cells_per_donor  <- 10
min_donors_per_group <- 3

deg_fdr_cut <- 0.05  # ONLY FDR

# Hallmark GSEA thresholds
hall_minGS <- 15
hall_maxGS <- 500

# GO GSEA thresholds
go_ont   <- "BP"
go_minGS <- 10
go_maxGS <- 800
go_pcut  <- 0.05

# =======================
# Two explicit gene filters (replace filterByExpr)
# 1) donor-level coverage:
#    detected if CPM > 1
#    keep if in SAT >= max(10, 20%) AND in VAT >= max(10, 20%)
# 2) cell-level detection:
#    keep if >= 2% nuclei have counts>0
# =======================
donor_cpm_cutoff <- 1
min_cell_pct     <- 0.01
min_donor_abs    <- 5
min_donor_frac   <- 0.20

# =======================
# Paths
# =======================
metadata_path <- "../../../../../../../data/08.integration_snRNA_object/01.freeze_version/metadata_after_libraryQC_after_ctQC_freeze.tsv"

#base_count_dir <- "~/scratch/project/00.adipose/data/07.full_scale_snRNA/08.integration/all/scvi/03.count_matrix"
base_count_dir <- "../../../../../../../data/08.integration_snRNA_object/scvi/03.count_matrix"
count_path_sat <- file.path(base_count_dir, "SAT", "01.raw_count", paste0(celltype, ".raw_count.tsv"))
count_path_vat <- file.path(base_count_dir, "VAT", "01.raw_count", paste0(celltype, ".raw_count.tsv"))

save_root <- file.path(
  "../../../../../../../data/21.differential_expressed_genes/01.pseudobulk/04.only_VAT_vs_SAT_pseudobulk/",
  paste0("BMIgrp_", bmi_group)
)

deg_df_dir    <- file.path(save_root, "01.DEG", "01.table")
deg_plot_dir  <- file.path(save_root, "01.DEG", "02.plot", celltype)
gsea_df_dir   <- file.path(save_root, "02.GSEA", "01.table")
gsea_plot_dir <- file.path(save_root, "02.GSEA", "02.plot", celltype)

dir.create(deg_df_dir,    recursive = TRUE, showWarnings = FALSE)
dir.create(deg_plot_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(gsea_df_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(gsea_plot_dir, recursive = TRUE, showWarnings = FALSE)

box_dir <- file.path(deg_plot_dir, "boxplots")
dir.create(box_dir, showWarnings = FALSE)

feature_dir <- file.path(deg_plot_dir, "featureplots_top20")
dir.create(feature_dir, showWarnings = FALSE)

# =======================
# Sanity check paths
# =======================
message("[Args] celltype=", celltype, " | bmi_group=", bmi_group)
message("[Path] SAT count = ", count_path_sat)
message("[Path] VAT count = ", count_path_vat)

if (!file.exists(count_path_sat)) stop("SAT count_path does not exist: ", count_path_sat)
if (!file.exists(count_path_vat)) stop("VAT count_path does not exist: ", count_path_vat)

# =======================
# Load metadata
# =======================
metadata <- fread(metadata_path, sep = "\t", header = TRUE, data.table = FALSE)
metadata$BMI_group <- ifelse(metadata$BMI < 27.5,"non_obese","obese")
# =======================
# Protein-coding gene set from GTF-derived table
# =======================
anno_file <- "../../../../../../../../00.reference/hg38/cellranger_reference/gene_biotype_from_gtf.tsv"
anno <- data.table::fread(anno_file)

pc_genes <- anno[gene_type == "protein_coding" & !is.na(gene_name) & gene_name != "", gene_name]
pc_genes <- unique(pc_genes)

# =======================
# Helper: read raw_count.tsv -> Seurat
# =======================
read_count_as_seurat <- function(fp) {
  m <- fread(fp, sep = "\t", header = TRUE, check.names = FALSE, data.table = FALSE)
  first_col <- colnames(m)[1]
  rownames(m) <- m[[first_col]]
  m[[first_col]] <- NULL
  m <- t(as.matrix(m))  # genes x cells
  CreateSeuratObject(counts = m)
}

# =======================
# Helper: Hallmark GSEA
# =======================
run_hallmark_gsea <- function(res_df, out_prefix, out_table_dir, out_plot_dir) {

  gene_list <- sign(res_df$logFC) * (-log10(res_df$PValue))
  names(gene_list) <- res_df$gene
  gene_list <- gene_list[!is.na(gene_list) & is.finite(gene_list)]
  gene_list <- gene_list[order(abs(gene_list), decreasing = TRUE)]
  gene_list <- gene_list[!duplicated(names(gene_list))]
  gene_list <- sort(gene_list, decreasing = TRUE)

  hallmark_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>%
    dplyr::select(gs_name, gene_symbol) %>%
    dplyr::distinct() %>%
    dplyr::rename(term = gs_name, gene = gene_symbol)

  eh <- GSEA(
    geneList      = gene_list,
    TERM2GENE     = hallmark_t2g,
    pAdjustMethod = "fdr",
    pvalueCutoff  = 1,
    minGSSize     = hall_minGS,
    maxGSSize     = hall_maxGS,
    verbose       = FALSE,
    seed          = TRUE,
    eps           = 0
  )

  if (is.null(eh) || nrow(as.data.frame(eh@result)) == 0) {
    message("[!] No Hallmark terms returned: ", out_prefix)
    return(invisible(NULL))
  }

  res_h <- as.data.frame(eh@result) %>%
    mutate(
      Direction = ifelse(NES > 0, "Up in VAT", "Up in SAT"),
      Description = gsub("^HALLMARK_", "", Description),
      Description = gsub("_", " ", Description)
    )

  write_csv(res_h, file.path(out_table_dir, paste0(out_prefix, ".Hallmark_GSEA.all.csv")))
  sig_h <- res_h %>% filter(p.adjust < 0.05)
  write_csv(sig_h, file.path(out_table_dir, paste0(out_prefix, ".Hallmark_GSEA.sig_FDR0.05.csv")))

  plot_src <- if (nrow(sig_h) > 0) sig_h else res_h
  topn <- 15
  plot_df <- plot_src %>%
    group_by(Direction) %>%
    arrange(p.adjust, .by_group = TRUE) %>%
    slice_head(n = topn) %>%
    ungroup() %>%
    mutate(Description = forcats::fct_reorder(Description, NES))

  if (nrow(plot_df) > 0) {
    p <- ggplot(plot_df, aes(x = NES, y = Description, fill = p.adjust)) +
      geom_col(width = 0.7) +
      scale_fill_gradient(high = "#377eb8", low = "#e41a1c",
                          name = "FDR (p.adjust)",
                          guide = guide_colorbar(reverse = TRUE)) +
      labs(
        x = "NES",
        y = NULL,
        title = paste0(celltype, " | Hallmark GSEA (VAT vs SAT) | ", bmi_group),
        subtitle = if (nrow(sig_h) > 0) "Terms shown: FDR<0.05" else "No FDR<0.05 terms; showing top by FDR"
      ) +
      theme_minimal(base_size = 11) +
      theme(plot.title = element_text(size = 14, face = "bold"))

    pdf(file.path(out_plot_dir, paste0(out_prefix, ".top15.Hallmark_GSEA.pdf")), width = 10, height = 8)
    print(p)
    dev.off()
  }

  invisible(TRUE)
}

# =======================
# Main
# =======================
deg_analysis <- function(celltype) {

  message("[] Read SAT + VAT and merge: ", celltype)

  seu_sat <- read_count_as_seurat(count_path_sat)
  seu_vat <- read_count_as_seurat(count_path_vat)
  seu <- merge(seu_sat, y = seu_vat)

  # ---- add metadata ----
  meta_sub <- metadata[metadata$barcode %in% rownames(seu@meta.data), ]
  rownames(meta_sub) <- meta_sub$barcode
  seu <- AddMetaData(seu, meta_sub)
  seu <- JoinLayers(seu)
  # ---- find BMI group column (auto) ----
  bmi_cols_try <- c("BMI_group")
  bmi_col <- bmi_cols_try[bmi_cols_try %in% colnames(seu@meta.data)]
  if (length(bmi_col) == 0) {
    stop(
      "Cannot find BMI-group column in metadata. Tried: ",
      paste(bmi_cols_try, collapse = ", "),
      "\nAvailable columns include:\n",
      paste(head(colnames(seu@meta.data), 80), collapse = ", ")
    )
  }
  bmi_col <- bmi_col[1]
  message("[BMI] using column: ", bmi_col)

  # ---- subset Asians + Depot + BMI group ----
  seu <- subset(
    seu,
    Race %in% c("Malay", "Chinese", "Indian") &
      Depot %in% c("SAT", "VAT") &
      BMI_group == bmi_group
  )
  if (ncol(seu) == 0) {
    message("[!] Skip: no cells after filtering Asians + Depot + BMI group=", bmi_group)
    return(invisible(NULL))
  }

  seu$Depot <- factor(seu$Depot, levels = c("SAT", "VAT"))

  # =======================
  # Donor filtering: each depot separately, then intersect
  # =======================
  df_cnt <- seu@meta.data %>%
    count(donor_id, Depot, name = "n")

  donors_sat_ok <- df_cnt %>% filter(Depot == "SAT", n >= min_cells_per_donor) %>% pull(donor_id)
  donors_vat_ok <- df_cnt %>% filter(Depot == "VAT", n >= min_cells_per_donor) %>% pull(donor_id)
  donors_keep <- intersect(donors_sat_ok, donors_vat_ok)

  message("[DonorFilter] SAT donors ok: ", length(unique(donors_sat_ok)))
  message("[DonorFilter] VAT donors ok: ", length(unique(donors_vat_ok)))
  message("[DonorFilter] donors kept (intersect): ", length(unique(donors_keep)))

  if (length(donors_keep) < min_donors_per_group) {
    message("[!] Skip: insufficient paired donors after filtering. kept=", length(donors_keep))
    return(invisible(NULL))
  }

  seu <- subset(seu, donor_id %in% donors_keep)

  # =======================
  # Pseudobulk per donor x Depot
  # =======================
  pb_counts <- AggregateExpression(
    seu,
    assays = "RNA",
    group.by = c("donor_id", "Depot"),
    slot = "counts"
  )$RNA
  df_counts <- as.data.frame(pb_counts)

  # =======================
  # Pseudobulk metadata
  # =======================
  meta_donor <- seu[[]] %>%
    dplyr::select(donor_id, Age, BMI, Race, Gender) %>%
    dplyr::distinct(donor_id, .keep_all = TRUE)

  meta_cell <- seu[[]] %>%
    dplyr::select(donor_id, Depot, Age, BMI, Race, Gender) %>%
    mutate(Depot = factor(Depot, levels = c("SAT", "VAT")))

  meta_pb <- meta_cell %>% distinct(donor_id, Depot, .keep_all = TRUE)

  df_ncells <- seu@meta.data %>%
    count(donor_id, Depot, name = "n_cells")

  meta_pb <- meta_pb %>%
    left_join(df_ncells, by = c("donor_id", "Depot")) %>%
    mutate(
      n_cells = ifelse(is.na(n_cells), 0, n_cells),
      log_n_cells = log10(n_cells + 1)
    )

  meta_pb$pb_id <- paste(meta_pb$donor_id, meta_pb$Depot, sep = "_")
  pb_ids <- colnames(df_counts)
  meta_pb <- meta_pb[match(pb_ids, meta_pb$pb_id), ]
  stopifnot(all(meta_pb$pb_id == pb_ids))

  meta_pb$Gender   <- factor(meta_pb$Gender)
  meta_pb$Race     <- factor(meta_pb$Race)
  meta_pb$Depot    <- relevel(factor(meta_pb$Depot), ref = "SAT")
  meta_pb$donor_id <- factor(meta_pb$donor_id)

  # =======================
  # edgeR (paired VAT vs SAT)
  # =======================
  dge <- DGEList(counts = df_counts)
  dge <- calcNormFactors(dge)

  design <- model.matrix(~ donor_id + Depot, data = meta_pb)


  # =======================
  # Two explicit gene filters (NO filterByExpr)
  # =======================
  # Filter 1 (cell-level): >=1% nuclei have counts>0
  cnt_cell <- LayerData(seu, assay = "RNA", layer = "counts")
  pct_cells <- Matrix::rowMeans(cnt_cell > 0)
  keep_cell <- pct_cells >= min_cell_pct

  # Filter 2 (donor-level): CPM>1 detected in enough SAT and VAT pseudobulk samples
  cpm_mat <- cpm(dge, normalized.lib.sizes = TRUE, log = FALSE)
  write.table(cpm_mat,paste0("~/scratch/project/00.adipose/script/01.adipose_sqtl/21.differential_expressed_genes/by_comparisons/depot_difference_by_celltype_by_BMI_group/",bmi_group,"/",bmi_group,".",celltype,".tsv"),sep = "\t",quote= F)
  is_sat <- meta_pb$Depot == "SAT"
  is_vat <- meta_pb$Depot == "VAT"

  n_sat <- sum(is_sat); n_vat <- sum(is_vat)
  k_sat <- max(min_donor_abs, ceiling(min_donor_frac * n_sat))
  k_vat <- max(min_donor_abs, ceiling(min_donor_frac * n_vat))

  det_sat <- rowSums(cpm_mat[, is_sat, drop = FALSE] > donor_cpm_cutoff) >= k_sat
  det_vat <- rowSums(cpm_mat[, is_vat, drop = FALSE] > donor_cpm_cutoff) >= k_vat
  keep_donor <- det_sat | det_vat # union genes

  genes_pb <- rownames(dge$counts)
  genes_ok <- intersect(genes_pb, names(keep_cell))
  keep <- keep_donor[genes_ok] & keep_cell[genes_ok]
  keep[is.na(keep)] <- FALSE
  
 # Filter 3 (biotype-level): keep only protein_coding genes
 keep_pc <- genes_ok %in% pc_genes
 keep <- keep & keep_pc

 message("[GeneFilter] kept after protein_coding: ", sum(keep))

  message("[GeneFilter] genes before: ", nrow(dge$counts))
  message("[GeneFilter] kept: ", sum(keep),
        " | donor coverage: CPM>", donor_cpm_cutoff,
        " in SAT>=", k_sat, " OR VAT>=", k_vat,
        " | cell detection>=", min_cell_pct*100, "% nuclei",
        " | protein_coding only")

  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge)

  dge <- estimateDisp(dge, design, robust = TRUE)
  fit <- glmQLFit(dge, design, robust = TRUE)

  # VAT vs SAT
  qlf <- glmQLFTest(fit, coef = "DepotVAT")

  res <- topTags(qlf, n = Inf)$table
  res$gene <- rownames(res)

  # save full & sig (ONLY FDR)
  write.table(
    res,
    file.path(deg_df_dir, paste0(celltype, ".", bmi_group, ".paired.edgeR.all.tsv")),
    sep = "\t", quote = FALSE, row.names = TRUE
  )

  deg_sig <- res %>% filter(FDR < deg_fdr_cut)
  write.table(
    deg_sig,
    file.path(deg_df_dir, paste0(celltype, ".", bmi_group, ".paired.edgeR.sig_FDR", deg_fdr_cut, ".tsv")),
    sep = "\t", quote = FALSE, row.names = TRUE
  )

  message("[] edgeR done: ", celltype, " | BMI group=", bmi_group,
          " | sig genes (FDR<", deg_fdr_cut, ") = ", nrow(deg_sig))
   # =======================
  # Per-gene boxplot (sig genes)
  # =======================
  if (nrow(deg_sig) > 0) {
    norm_cpm <- cpm(dge, normalized.lib.sizes = TRUE, log = FALSE)
    log_norm_cpm <- log2(norm_cpm + 1)

    plot_data <- reshape2::melt(log_norm_cpm)
    colnames(plot_data) <- c("Gene", "donor_id", "Expression")
    head(plot_data)
    cols <- c("SAT" = "#BBCFC3", "VAT" = "#FBCE9E")
    ylab_txt <- "log2(CPM+1)"
    for (g in deg_sig$gene) {
      expr_tmp <- plot_data %>%
  filter(Gene == g) %>%
  mutate(
    Depot    = sub("^.*_(SAT|VAT)$", "\\1", donor_id),
    donor_id = sub("_(SAT|VAT)$", "", donor_id),
    Depot    = factor(Depot, levels = c("SAT", "VAT"))
  )
      st <- deg_sig[g, , drop = FALSE]
      fdr_txt <- formatC(st$FDR[1], format = "e", digits = 2)
      lfc_txt <- sprintf("%.2f", st$logFC[1])
      subtitle_txt <- paste0("logFC = ", lfc_txt, "   FDR = ", fdr_txt)
    p <- ggplot(expr_tmp, aes(x = Depot, y = Expression)) +
        geom_boxplot(aes(fill = Depot), width = 0.5, alpha = 0.7, outlier.shape = NA) +
        geom_line(aes(group = donor_id), color = "grey55", linewidth = 0.25, alpha = 0.8) +
        geom_point(aes(color = Depot), size = 1.8, alpha = 0.9) +
        scale_fill_manual(values = cols) +
        scale_color_manual(values = cols) +
        labs(
          title = paste0(celltype, ": ", g),
          subtitle = subtitle_txt,
          x = "Depot",
          y = ylab_txt
        ) +
        theme_bw(base_size = 13) +
        theme(legend.position = "none", plot.subtitle = element_text(size = 6))

      pdf(file.path(box_dir, paste0(g, ".boxplot.pdf")), width = 3.2, height = 3.2)
      print(p)
      dev.off()
    }
  }
  # =======================
  # Volcano (ONLY FDR)
  # =======================
  res2 <- res %>%
    mutate(group = case_when(
      FDR < deg_fdr_cut & logFC > 1  ~ "Up in VAT",
      FDR < deg_fdr_cut & logFC < -1  ~ "Up in SAT",
      TRUE                           ~ "NotSig"
    ))

  top_up   <- res2 %>% filter(group == "Up in VAT") %>% arrange(FDR) %>% slice_head(n = 15)
  top_down <- res2 %>% filter(group == "Up in SAT") %>% arrange(FDR) %>% slice_head(n = 15)
  top_label <- bind_rows(top_up, top_down)

  pdf(file.path(deg_plot_dir, paste0(celltype, ".", bmi_group, ".volcano.VAT_vs_SAT.FDR", deg_fdr_cut, ".pdf")),
      width = 6, height = 4)
  print(
    ggplot(res2, aes(x = logFC, y = -log10(FDR))) +
      geom_point(data = subset(res2, group == "NotSig"), color = "grey80", size = 1, alpha = 0.6) +
      geom_point(data = subset(res2, group == "Up in VAT"), color = "red", size = 1, alpha = 0.7) +
      geom_point(data = subset(res2, group == "Up in SAT"), color = "blue", size = 1, alpha = 0.7) +
      geom_text_repel(data = top_label, aes(label = gene), size = 3,
                      max.overlaps = Inf, box.padding = 0.4, point.padding = 0.3) +
      geom_hline(yintercept = -log10(deg_fdr_cut), linetype = "dashed") +
      geom_vline(xintercept = 0, linetype = "solid") +
      theme_bw() +
      labs(
        title = paste0(celltype, " | VAT vs SAT (paired) | ", bmi_group),
        x = "log2 Fold Change (VAT / SAT)",
        y = "-log10(FDR)"
      )
  )
  dev.off()

  # =======================
  # FeaturePlot: top 20 genes (each gene one PDF, 6x4)
  # =======================
  set.seed(123)
  plot_genes <- deg_sig %>% arrange(FDR) %>% slice_head(n = 20) %>% pull(gene)
  if (length(plot_genes) < 1) plot_genes <- res %>% arrange(FDR) %>% slice_head(n = 20) %>% pull(gene)

  message("[FeaturePlot] Running quick PCA/UMAP for plotting...")
  DefaultAssay(seu) <- "RNA"
  seu <- NormalizeData(seu, verbose = FALSE)
  seu <- FindVariableFeatures(seu, nfeatures = 2000, verbose = FALSE)
  seu <- ScaleData(seu, features = VariableFeatures(seu), verbose = FALSE)
  seu <- RunPCA(seu, features = VariableFeatures(seu), npcs = 30, verbose = FALSE)
  seu <- RunUMAP(seu, dims = 1:30, verbose = FALSE)

  for (g in plot_genes) {
    if (!g %in% rownames(seu)) next
    pdf(file.path(feature_dir, paste0(g, ".FeaturePlot.", bmi_group, ".pdf")), width = 6, height = 4)
    print(
      FeaturePlot(seu, features = g, reduction = "umap", cols = c("lightgrey", "red")) +
        ggtitle(paste0(celltype, " | ", g, " | ", bmi_group))
    )
    dev.off()
  }

  # =======================
  # GO GSEA (table)
  # =======================
  gene_list_go <- sign(res$logFC) * (-log10(res$PValue))
  names(gene_list_go) <- res$gene
  gene_list_go <- gene_list_go[!is.na(gene_list_go) & is.finite(gene_list_go)]
  gene_list_go <- gene_list_go[order(abs(gene_list_go), decreasing = TRUE)]
  gene_list_go <- gene_list_go[!duplicated(names(gene_list_go))]
  gene_list_go <- sort(gene_list_go, decreasing = TRUE)

  ego <- gseGO(
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
  )

  if (!is.null(ego) && nrow(as.data.frame(ego@result)) > 0) {
    res_go <- as.data.frame(ego@result)
    write_csv(res_go, file.path(gsea_df_dir, paste0(celltype, ".", bmi_group, "_gseGO_", go_ont, ".csv")))
  } else {
    message("[!] No significant gseGO terms.")
  }

  # =======================
  # Hallmark GSEA (all + sig + plot)
  # =======================
  run_hallmark_gsea(
    res_df = res,
    out_prefix = paste0(celltype, ".", bmi_group, ".VAT_vs_SAT"),
    out_table_dir = gsea_df_dir,
    out_plot_dir = gsea_plot_dir
  )

  message("[] ALL DONE: ", celltype, " | BMI group=", bmi_group)
  invisible(TRUE)
}

deg_analysis(celltype = celltype)

