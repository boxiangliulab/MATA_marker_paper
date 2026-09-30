## GTEx replication of MAATA eQTLs
## Assumes the following files already exist:
##   All_eQTL_SAT_and_VAT.tsv            (maata_anno)
##   All_eQTL_SAT_and_VAT_with_beta.tsv  (merged, contains MAATA BETA per lead SNP)

library(arrow)
library(data.table)
library(ggplot2)
library(ggrastr)

colors <- c("Adipocyte" = "#56BF95",
            "Early_preadipocyte"= "#c7e9c0",
            "Committed_preadipocyte" = "#a1d99b",
            "Areg" = "#CCCC4D",
            "IGFBP2+_cell" = "#74c476",
            "Mesothelial" ="#DAE220",
            "AEC" ="#F0B6DA",
            "CEC" = "#DE77AE",
            "VEC" = "#C51B7D",
            "LEC" = "#8E0052",
            "Pericyte"="#fcc5c0",
            "SMC2"="#fa9fb5",
            "SMC"="#f768a1",
            "SMC1"="#f768a1",
            "Endometrial"="#883186",
            "cDC1" = "#a6bddb",
            "cDC2" = "#74a9cf",
            "Activate_DC" = "#3690c0",
            "pDC"="#0570b0",
            "Classical_monocyte" = "#243494",
            "Non_classical_monocyte"="#061C58",
            "Intermediate_monocyte" = "#1C1894",
            "TIM4+_ATM"="#bcbddc",
            "PVM"="#9e9ac8",
            "LAM" = "#807dba",
            "TIM4+_CD11c+_ATM" = "#6a51a3",
            "B" = "#fdae61",
            "Plasmablast"="#993404",
            "Plasma_cell"="#f46d43",
            "Treg" = "#A7D5F7",
            "CD4+_T"="#9BB9FB",
            "CD8+_T"="#5087F3",
            "CD56dim_CD16+_NK" = "#4157E4",
            "CD56+_NK" = "#253CD1",
            "Mast"="#44B3E1")

# ── paths ──────────────────────────────────────────────────────────────────────
# Project root is 8 levels up from this script
args0        <- commandArgs(trailingOnly = FALSE)
script_path  <- sub("--file=", "", args0[grep("--file=", args0)])
project_root <- normalizePath(file.path(dirname(normalizePath(script_path)),
                                        paste(rep("..", 8), collapse = "/")))

maata_base <- file.path(project_root, "data/09.eQTL/03.Saige_QTL/02.output/03.formal_calling/ancestry_specific_version")
out_dir    <- file.path(project_root, "data/09.eQTL/03.Saige_QTL/02.output/04.replication")
gtex_base  <- "/path/to/GTEx_Analysis_v10_eQTL_updated"  # GTEx v10; downloaded from https://gtexportal.org
gtf_path   <- normalizePath(file.path(project_root, "../00.reference/hg38/cellranger_reference/refdata-gex-GRCh38-2020-A/genes/genes.gtf"), mustWork = FALSE)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── load pre-computed MAATA tables ─────────────────────────────────────────────
maata_anno <- fread(file.path(maata_base, "All_eQTL_SAT_and_VAT.tsv"))
merged     <- fread(file.path(maata_base, "All_eQTL_SAT_and_VAT_with_beta.tsv"))
setDT(maata_anno)
setDT(merged)

# ── GTEx gene_id -> gene_name mapping from GTF ────────────────────────────────
gtf_dt <- fread(gtf_path, sep = "\t", header = FALSE, quote = "",
                fill = TRUE, comment.char = "#", data.table = TRUE)
stopifnot(ncol(gtf_dt) >= 9)
setnames(gtf_dt, paste0("V", seq_len(ncol(gtf_dt))))
gtf_gene <- gtf_dt[V3 == "gene", .(attr = V9)]
gtf_gene[, gene_id   := sub('.*gene_id "([^"]+)".*',   '\\1', attr)]
gtf_gene[, gene_name := sub('.*gene_name "([^"]+)".*', '\\1', attr)]
anno <- unique(gtf_gene[, .(gene_id, gene_name)])
setDT(anno); setkey(anno, gene_id)

# ── GTEx replication function ─────────────────────────────────────────────────
gtex_rep <- function(gtex, depot) {

  gtex$variant_id <- gsub("_", ":", sub("_b38$", "", gtex$variant_id))
  setDT(gtex)
  gtex[, gene_id_nov := sub("\\.[0-9]+$", "", gene_id)]

  gtex2 <- merge(gtex, anno, by.x = "gene_id_nov", by.y = "gene_id", all.x = TRUE)
  gtex2[, eQTL := paste0(gene_name, "-", variant_id)]

  maata_tmp <- maata_anno[depot == depot]

  # (A) depot-level replication rate
  eqtls_all <- unique(maata_tmp$eQTL)
  total_n_maata <- length(eqtls_all)
  total_n_maata_replicate <- sum(eqtls_all %in% gtex2$eQTL)

  summary_depot <- data.frame(
    source = "GTEx",
    depot  = depot,
    total_n_maata_eqtl           = total_n_maata,
    total_n_maata_replicate_eqtl = total_n_maata_replicate,
    replicate_ratio = if (total_n_maata == 0) NA_real_ else total_n_maata_replicate / total_n_maata
  )

  # (B) celltype-level overlap
  cts <- sort(unique(maata_tmp$celltype))
  out <- vector("list", length(cts))
  for (k in seq_along(cts)) {
    ct <- cts[k]
    eqtls_ct  <- unique(maata_tmp[celltype == ct, eQTL])
    n_total_ct <- length(eqtls_ct)
    n_overlap  <- sum(eqtls_ct %in% gtex2$eQTL)
    out[[k]] <- data.frame(
      source    = "GTEx",
      depot     = depot,
      celltype  = ct,
      total_maata_eqtl_celltype = n_total_ct,
      overlap_eqtl_celltype     = n_overlap,
      overlap_ratio_celltype    = if (n_total_ct == 0) NA_real_ else n_overlap / n_total_ct
    )
  }

  list(summary_depot   = summary_depot,
       overlap_celltype = rbindlist(out))
}

# SAT
gtex_sat <- read_parquet(file.path(gtex_base, "Adipose_Subcutaneous.v10.eQTLs.signif_pairs.parquet"))
res_sat   <- gtex_rep(gtex_sat, "SAT")

# VAT
gtex_vat <- read_parquet(file.path(gtex_base, "Adipose_Visceral_Omentum.v10.eQTLs.signif_pairs.parquet"))
res_vat   <- gtex_rep(gtex_vat, "VAT")

# save replication summary tables
rep_summary <- rbind(res_sat$summary_depot,   res_vat$summary_depot)
ct_overlap  <- rbind(res_sat$overlap_celltype, res_vat$overlap_celltype)
fwrite(rep_summary, file.path(out_dir, "GTEx_replication_depot_summary.tsv"),   sep = "\t")
fwrite(ct_overlap,  file.path(out_dir, "GTEx_replication_celltype_overlap.tsv"), sep = "\t")

# ── build rep_dt: MAATA BETA vs GTEx slope ────────────────────────────────────
# preprocess GTEx (both depots together)
gtex_sat$depot <- "SAT"; gtex_vat$depot <- "VAT"
gtex_all <- rbind(gtex_sat, gtex_vat)
setDT(gtex_all)
gtex_all[, variant_id_clean := gsub("_", ":", sub("_b38$", "", variant_id))]
gtex_all[, gene_id_nov      := sub("\\.[0-9]+$", "", gene_id)]
gtex_all2 <- merge(gtex_all, anno, by.x = "gene_id_nov", by.y = "gene_id", all.x = TRUE)
gtex_all2[, eQTL := paste0(gene_name, "-", variant_id_clean)]

# merged has eqtl_key (= eQTL from maata_anno) and BETA (our effect size)
maata_beta <- merged[!is.na(BETA), .(depot, celltype, eQTL = eqtl_key, our_BETA = BETA)]

# GTEx slope lookup (one row per depot+eQTL; take first if duplicates)
gtex_slope_dt <- unique(gtex_all2[!is.na(slope), .(depot, eQTL, gtex_slope = slope)])

rep_dt <- merge(maata_beta, gtex_slope_dt, by = c("depot", "eQTL"))
setDT(rep_dt)

# ── beta correlation scatter plots ────────────────────────────────────────────
min_repli <- 50

keep_ct <- rep_dt[!is.na(our_BETA) & !is.na(gtex_slope),
                  .N, by = .(depot, celltype)][N >= min_repli]
rep_dt_f <- rep_dt[keep_ct, on = .(depot, celltype)]
rep_dt_f[, celltype := factor(celltype, levels = names(colors))]

plot_one_depot <- function(dep) {
  df <- rep_dt_f[depot == dep & !is.na(our_BETA) & !is.na(gtex_slope)]
  if (nrow(df) == 0) return(invisible(NULL))
  df[, celltype := droplevels(celltype)]

  p <- ggplot(df, aes(x = our_BETA, y = gtex_slope)) +
    geom_point_rast(aes(color = celltype), size = 0.8, alpha = 0.7, dpi = 600) +
    geom_hline(yintercept = 0, linewidth = 0.25) +
    geom_vline(xintercept = 0, linewidth = 0.25) +
    scale_color_manual(values = colors, drop = TRUE, guide = "none") +
    facet_wrap(~ celltype, scales = "free", ncol = 6) +
    theme_bw(base_size = 10) +
    theme(
      strip.text   = element_text(size = 8),
      plot.title   = element_text(face = "bold"),
      plot.caption = element_text(size = 7, hjust = 0, margin = margin(t = 8))
    ) +
    labs(
      title = paste0("MAATA BETA vs GTEx slope (", dep, "), celltypes with n \u2265 ", min_repli),
      x = "MAATA eQTL effect (BETA)",
      y = "GTEx eQTL effect (slope)",
      caption = paste0(
        "Figure. Scatter plots showing the correlation of eQTL effect sizes (BETA) between MAATA snRNA-seq ",
        "and GTEx v10 ", dep, " adipose tissue (slope). Each point represents a shared significant eQTL ",
        "(gene\u2013variant pair) detected in both datasets. Only cell types with \u2265 ", min_repli,
        " overlapping eQTLs are shown. Points are coloured by cell type. ",
        "Horizontal and vertical dashed lines indicate zero effect."
      )
    )

  ggsave(file.path(out_dir, paste0("scatter_our_vs_gtex_", dep, "_facet_celltypes.min", min_repli, ".pdf")),
         p, width = 12, height = 6, units = "in", limitsize = FALSE)

  fwrite(df, file.path(out_dir, paste0("replication_intersection_", dep, ".min", min_repli, ".tsv")), sep = "\t")
}

for (dep in sort(unique(rep_dt_f$depot))) plot_one_depot(dep)
