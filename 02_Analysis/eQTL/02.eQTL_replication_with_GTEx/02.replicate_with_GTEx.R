## AdipoExpress is converted into hg38
## GTEx is hg38 originally
## MAATA eQTL: lead SNP - eGene - cell type - depot

# ── project root (8 levels up from this script) ────────────────────────────────
args0        <- commandArgs(trailingOnly = FALSE)
script_path  <- sub("--file=", "", args0[grep("--file=", args0)])
project_root <- normalizePath(file.path(dirname(normalizePath(script_path)),
                                        paste(rep("..", 8), collapse = "/")))
# GTEx v10 eQTL results; downloaded from https://gtexportal.org
gtex_base    <- "/path/to/GTEx_Analysis_v10_eQTL_updated"

library(arrow)
library(data.table)
library(qvalue)
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
## Save the association files of MAATA first 
n_overlap = 50

anno_file <- normalizePath(file.path(project_root, "../00.reference/hg38/cellranger_reference/gene_biotype_from_gtf.tsv"), mustWork = FALSE)
anno <- data.table::fread(anno_file)
setDT(anno)
## Merge of all celltype - depot gene-SNP pair (only egene)
# SAT
celltypes <- read.table(file.path(project_root, "script/09.eQTL/03.Saige_QTL/04.QTL_calling/02.output/03.formal_calling/ancestry_specific_version/SAT/00.celltype_to_test.list"))
base_dir <- file.path(project_root, "data/09.eQTL/03.Saige_QTL/02.output/03.formal_calling/ancestry_specific_version/SAT/01.egene_list")

lst <- list()
k <- 0

for (ct in celltypes$V1) {
  f <- file.path(base_dir, paste0("SAT.", ct, ".egene_list.tsv"))

  if (!file.exists(f)) {
    message("[SKIP missing] ", f)
    next
  }

  tmp <- tryCatch(
    read.table(f, sep = "\t", header = TRUE, stringsAsFactors = FALSE),
    error = function(e) {
      message("[SKIP read error] ", f, " :: ", e$message)
      return(NULL)
    }
  )

  if (is.null(tmp) || nrow(tmp) == 0) {
    message("[SKIP empty] ", f)
    next
  }

  tmp$celltype <- ct  
  k <- k + 1
  lst[[k]] <- tmp
}
sat_df <- if (length(lst) == 0) data.frame() else do.call(rbind, lst)
# VAT
celltypes <- read.table(file.path(project_root, "script/09.eQTL/03.Saige_QTL/04.QTL_calling/02.output/03.formal_calling/ancestry_specific_version/VAT/00.celltype_to_test.list"))
base_dir <- file.path(project_root, "data/09.eQTL/03.Saige_QTL/02.output/03.formal_calling/ancestry_specific_version/VAT/01.egene_list")

lst <- list()
k <- 0

for (ct in celltypes$V1) {
  f <- file.path(base_dir, paste0("VAT.", ct, ".egene_list.tsv"))

  if (!file.exists(f)) {
    message("[SKIP missing] ", f)
    next
  }

  tmp <- tryCatch(
    read.table(f, sep = "\t", header = TRUE, stringsAsFactors = FALSE),
    error = function(e) {
      message("[SKIP read error] ", f, " :: ", e$message)
      return(NULL)
    }
  )

  if (is.null(tmp) || nrow(tmp) == 0) {
    message("[SKIP empty] ", f)
    next
  }

  tmp$celltype <- ct 
  k <- k + 1
  lst[[k]] <- tmp
}

vat_df <- if (length(lst) == 0) data.frame() else do.call(rbind, lst)
maata_eqtl <- rbind(vat_df,sat_df)
maata_eqtl$gene <- gsub("_","-",maata_eqtl$gene)
maata_eqtl$eQTL <- paste0(maata_eqtl$gene,"-",maata_eqtl$top_MarkerID)
setDT(maata_eqtl)

# add gene annotation in the data frame
anno_u <- unique(anno, by = "gene_name")
maata_anno <- merge(
  maata_eqtl,
  anno_u,
  by.x = "gene",
  by.y = "gene_name",
  all.x = TRUE
)

write.table(maata_anno, file.path(project_root, "data/09.eQTL/03.Saige_QTL/02.output/03.formal_calling/ancestry_specific_version/All_eQTL_SAT_and_VAT.tsv"), sep = "\t", quote=F, row.names=F)

## Save the association files of MAATA first

# Read sum-stats for one celltype/depot, keep only lead-SNP rows, return BETA + eQTL2 key
sum_stats_base <- file.path(project_root, "data/09.eQTL/03.Saige_QTL/02.output/03.formal_calling/ancestry_specific_version")

generate_assoc <- function(ct, dep) {
  f <- file.path(sum_stats_base, dep, "02.sum_stats",
                 paste0(ct, "_all.gene.cis.tsv"))
  if (!file.exists(f)) { message("[SKIP missing] ", f); return(data.table()) }
  dt <- tryCatch(fread(f, sep = "\t"),
                 error = function(e) { message("[SKIP] ", f, ": ", e$message); NULL })
  if (is.null(dt) || nrow(dt) == 0) return(data.table())

  # gene names in sum-stats use "_"; maata_anno uses "-"
  dt[, gene_dash := gsub("_", "-", gene)]
  dt[, eQTL2     := paste0(gene_dash, "-", MarkerID)]

  # keep only eGene lead SNPs (those already recorded in maata_anno)
  leads <- unique(maata_anno[depot == dep & celltype == ct, eQTL])
  dt_lead <- dt[eQTL2 %in% leads]

  dt_lead[, .(depot, celltype, gene = gene_dash, MarkerID, BETA, SE, p.value, eQTL2)]
}

# SAT
celltypes <- read.table(file.path(project_root, "script/09.eQTL/03.Saige_QTL/04.QTL_calling/02.output/03.formal_calling/ancestry_specific_version/SAT/00.celltype_to_test.list"))
sat_ass <- data.frame()
for (i in celltypes$V1){
        tmp <- generate_assoc(i,"SAT") 
        sat_ass <- rbind(sat_ass,tmp)
        }

# VAT
celltypes <- read.table(file.path(project_root, "script/09.eQTL/03.Saige_QTL/04.QTL_calling/02.output/03.formal_calling/ancestry_specific_version/VAT/00.celltype_to_test.list"))
vat_ass <- data.frame()
for (i in celltypes$V1){
        tmp <- generate_assoc(i,"VAT")
        vat_ass <- rbind(vat_ass,tmp)
}
ass_all <- rbind(vat_ass,sat_ass)

# merge
maata_anno[, eqtl_key := trimws(as.character(eQTL))]
ass_all[,   eqtl_key := trimws(as.character(eQTL2))]
merged <- ass_all[maata_anno, on=.(depot, celltype, eqtl_key)]
write.table(merged, file.path(project_root, "data/09.eQTL/03.Saige_QTL/02.output/03.formal_calling/ancestry_specific_version/All_eQTL_SAT_and_VAT_with_beta.tsv"), sep = "\t", quote=F, row.names=F)
########## GTEx ########
gtex <- read_parquet(file.path(gtex_base, "Adipose_Subcutaneous.v10.eQTLs.signif_pairs.parquet"))
gtf <- normalizePath(file.path(project_root, "../00.reference/hg38/cellranger_reference/refdata-gex-GRCh38-2020-A/genes/genes.gtf"), mustWork = FALSE)

gtf_dt <- fread(
  gtf,
  sep = "\t",
  header = FALSE,
  quote = "",
  fill = TRUE,
  comment.char = "#",   
  data.table = TRUE
)

stopifnot(ncol(gtf_dt) >= 9)

setnames(gtf_dt, paste0("V", seq_len(ncol(gtf_dt))))  

gtf_gene <- gtf_dt[V3 == "gene", .(attr = V9)]

gtf_gene[, gene_id := sub('.*gene_id "([^"]+)".*', '\\1', attr)]
gtf_gene[, gene_name := sub('.*gene_name "([^"]+)".*', '\\1', attr)]
anno <- unique(gtf_gene[, .(gene_id, gene_name)])
head(anno)

gtex_rep <- function(gtex, depot){

  # ---- GTEx preprocessing ----
  gtex$variant_id <- gsub("_", ":", sub("_b38$", "", gtex$variant_id))
  setDT(gtex)
  gtex[, gene_id_nov := sub("\\.[0-9]+$", "", gene_id)]

  setDT(anno)
  setkey(anno, gene_id)

  gtex2 <- merge(
    gtex, anno,
    by.x = "gene_id_nov", by.y = "gene_id",
    all.x = TRUE
  )

  gtex2[, eQTL := paste0(gene_name, "-", variant_id)]

  # ---- subset MAATA ----
  setDT(maata_anno)
  maata_tmp <- maata_anno[depot == depot]

  # ============================
  # (A) Depot-level replication rate
  # ============================
  eqtls_all <- unique(maata_tmp$eQTL)
  total_n_maata <- length(eqtls_all)
  total_n_maata_replicate <- sum(eqtls_all %in% gtex2$eQTL)

  summary_depot <- data.frame(
    source = "GTEx",
    depot = depot,
    total_n_maata_eqtl = total_n_maata,
    total_n_maata_replicate_eqtl = total_n_maata_replicate,
    replicate_ratio = if (total_n_maata == 0) NA_real_ else total_n_maata_replicate / total_n_maata
  )

  # ============================
  # (B) Celltype-level pi1
  # ============================
  celltypes <- sort(unique(maata_tmp$celltype))
  out <- vector("list", length(celltypes))

  for (k in seq_along(celltypes)) {
    ct <- celltypes[k]

   eqtls_ct <- unique(maata_tmp[celltype == ct, eQTL])
    n_total_ct <- length(eqtls_ct)

   n_overlap <- sum(eqtls_ct %in% gtex2$eQTL)
   # p <- gtex2[eQTL %in% eqtls_ct, pval_nominal]
   # p <- p[is.finite(p) & !is.na(p) & p > 0 & p <= 1]

   # pi1 <- NA_real_
   # if (n_overlap >= 50) {
   #   pi1 <- tryCatch(1 - pi0est(p)$pi0, error = function(e) NA_real_)
   # }

    out[[k]] <- data.frame(
      source = "GTEx",
      depot = depot,
      celltype = ct,
      total_maata_eqtl_celltype = n_total_ct,
      overlap_eqtl_celltype = n_overlap,
      overlap_ratio_celltype = if (n_total_ct == 0) NA_real_ else n_overlap / n_total_ct,
    #  pi1 = pi1
    )
  }

  pi1_by_celltype <- rbindlist(out)

  return(list(
    summary_depot = summary_depot,
    pi1_by_celltype = pi1_by_celltype
  ))
}

# SAT
gtex <- read_parquet(file.path(gtex_base, "Adipose_Subcutaneous.v10.eQTLs.signif_pairs.parquet"))
res_sat <- gtex_rep(gtex, "SAT")
sat_summary <- res_sat$summary_depot
#sat_pi1 <- res_sat$pi1_by_celltype

# VAT
gtex <- read_parquet(file.path(gtex_base, "Adipose_Visceral_Omentum.v10.eQTLs.signif_pairs.parquet"))
res_vat <- gtex_rep(gtex, "VAT")
vat_summary <- res_vat$summary_depot
#sat_pi1 <- res_vat$pi1_by_celltype

############ visualize the beta by cell type in each depot #########
sat_gtex <- read_parquet(file.path(gtex_base, "Adipose_Subcutaneous.v10.eQTLs.signif_pairs.parquet"))
sat_gtex$depot <- "SAT"
vat_gtex <- read_parquet(file.path(gtex_base, "Adipose_Visceral_Omentum.v10.eQTLs.signif_pairs.parquet"))
vat_gtex$depot="VAT"
gtex_all <- rbind(vat_gtex,sat_gtex)

out_dir <- file.path(project_root, "data/09.eQTL/03.Saige_QTL/02.output/04.replication")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

setDT(rep_dt)

min_repli <- 50

# filter celltypes with enough replicated points PER DEPOT
keep_ct <- rep_dt[!is.na(our_BETA) & !is.na(gtex_slope),
                  .N, by=.(depot, celltype)][N >= min_repli]

rep_dt_f <- rep_dt[keep_ct, on=.(depot, celltype)]

# (optional) make celltype a factor to enforce palette order
rep_dt_f[, celltype := factor(celltype, levels = names(colors))]

# plotting function (one figure per depot; facet all remaining celltypes)
plot_one_depot <- function(dep) {
  df <- rep_dt_f[depot == dep & !is.na(our_BETA) & !is.na(gtex_slope)]
  if (nrow(df) == 0) return(invisible(NULL))
# drop unused levels for clean faceting
  df[, celltype := droplevels(celltype)]

  p <- ggplot(df, aes(x = our_BETA, y = gtex_slope)) +
    geom_point_rast(aes(color = celltype), size = 0.8, alpha = 0.7, dpi = 600) +
    geom_hline(yintercept = 0, linewidth = 0.25) +
    geom_vline(xintercept = 0, linewidth = 0.25) +
    scale_color_manual(values = colors, drop = TRUE, guide = "none") +
    facet_wrap(~ celltype, scales = "free", ncol = 6) +
    theme_bw(base_size = 10) +
    theme(
      strip.text = element_text(size = 8),
      plot.title = element_text(face = "bold"),
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

  # also save filtered intersection table for this depot
  fwrite(df, file.path(out_dir, paste0("replication_intersection_", dep, ".min", min_repli, ".tsv")), sep = "\t")
}

# run for all depots present
for (dep in sort(unique(rep_dt_f$depot))) plot_one_depot(dep)




