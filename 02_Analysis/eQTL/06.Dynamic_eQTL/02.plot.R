library(data.table)
library(tidyverse)
library(edgeR)
# metadata
metadata = read.table("/project/01.adipose_sqtl/data/08.integration_snRNA_object/01.freeze_version/metadata_after_libraryQC_after_ctQC_freeze.tsv",sep = "\t", header = T)
metadata <- metadata[-which(metadata$Race == "Caucasian"),]
donor_meta <- unique(metadata[,c("barcode","donor_id")])

slingshot_cond_gene <- read.table("/slingshot/gene_list.txt")
inter_gene <- read.table("/project/01.adipose_sqtl/data/24.dynamic_qtl/01.dynamic_eQTL/01.input/SAT_VAT_intersect_egene.tsv")


# read genotype
geno_raw <- fread("/project/adipose/fullgenotype/Ancestry_specific/00.merge_genotype/genotype_adipose_after_imputation_296_merged_genome_rename_rsID_neur.tsv",sep="\t")
geno_long <- geno_raw %>% pivot_longer(cols = -1, names_to = "donor_id", values_to = "genotype")
geno_long <- data.table(geno_long) 
setDT(geno_long)
setkey(geno_long, SNP)

# generate significant dynamic eQTL
sig_dynamic_eqtl <- function(depot){
  input <- paste0("/project/adipose/dynamic_qtl/02.output/",depot,"/dynamic_eQTL/")
  filelist <- list.files(input, pattern = "\\.tsv$", full.names = TRUE)
  df <- do.call(rbind,lapply(filelist, function(f) {
    read.delim(f, stringsAsFactors = FALSE)
  })
  )
  df <- df %>%
    mutate(
      p_dynamic_fdr = p.adjust(p_dynamic, method = "BH")
    )
  return(df)
}

sat <- sig_dynamic_eqtl("SAT")
vat <- sig_dynamic_eqtl("VAT")
merge <- rbind(sat,vat)
merge$lead_eqtl <- paste0(merge$gene_id,"-",merge$snp_lead)
merge_sig <- merge[which(merge$p_dynamic_fdr <= 0.05),]
sig_genes <- unique(merge_sig$gene_id)
sig_genes
gene_anno <- fread("/project/adipose/dynamic_qtl/01.input/SAT_VAT_merged_eQTL.tsv") %>%
  select(gene, gene_type) %>%
  distinct()
merge_sig <- merge_sig %>%
  left_join(gene_anno, by = c("gene_id" = "gene"))

merge_sig <- merge_sig %>%
  filter(gene_type == "protein_coding")

merge_sig <- merge_sig %>%
  filter(snp_id == snp_lead)
top_df <- merge_sig %>%
  as.data.frame() %>%
  distinct(depot, gene_id, snp_lead, .keep_all = TRUE) %>%
  group_by(gene_id, snp_lead) %>%
  filter(n_distinct(depot) == 2) %>%
  ungroup()

top_df <- top_df %>%
  distinct(gene_id, snp_lead, .keep_all = TRUE)
#top_df <- merge_sig %>%
  #as.data.frame() %>%
 # distinct(depot, gene_id, snp_lead, .keep_all = TRUE) %>%
  #group_by(gene_id, snp_lead) %>%
  #filter(n() == 1) %>%   # 只保留只出现在一个 depot 的
  #ungroup()

# =========================================================
# 2. Take top100 unique gene-snp signals
# =========================================================

top100  <- top_df %>%
    as.data.frame() %>%
    arrange(p_dynamic_fdr, p_dynamic) %>%
    distinct(gene_id, snp_lead, .keep_all = TRUE) %>%
    slice_head(n = 60)
library(dplyr)

# 你给的强制加入列表
extra <- tribble(
  ~gene_id, ~snp_id,
  "DLGAP2","chr8:1141262:C:T",
  "DLGAP2","chr8:1152269:T:G",
  "DLGAP2","chr8:1162047:C:A",
  "PASK","chr2:241526411:C:A",
  "AMPH","chr7:38470929:A:T",
  "ZNF696","chr8:143275043:G:C",
  "AMPH","chr7:38460448:A:G",
  "CLTCL1","chr22:19276024:A:G",
  "PON3","chr7:95396288:G:A",
  "DIAPH1","chr5:141615674:C:T",
  "DIAPH1","chr5:141623077:C:A",
  "DIAPH1","chr5:141617123:T:C",
  "SAMD12","chr8:118619440:C:G",
  "ZNF334","chr20:46529811:G:C"
)

# 从 merge_sig 里抓这些行
extra_rows <- merge_sig %>%
  semi_join(extra, by = c("gene_id","snp_id"))

# 合并到 top100
top100_final <- bind_rows(top100, extra_rows) %>%
  distinct(gene_id, snp_id, depot, .keep_all = TRUE)

top100<-top100_final



# first gene tested
gene_list <- unique(top100$gene_id)

celltype <- c("Adipocyte","Committed_preadipocyte","Early_preadipocyte")
sub <- merge_sig[which(merge_sig$gene_id %in% gene_list),]
plot_pairs0 <- unique(sub[, c("gene_id", "snp_lead")])
plot_pairs0 <- as.data.frame(plot_pairs0)

plot_pairs <- merge(
  plot_pairs0,
  data.frame(depot = c("SAT", "VAT")),
  by = NULL
)

plot_pairs <- plot_pairs[, c("gene_id", "depot", "snp_lead")]

geno_sub <- geno_long[.(unique(plot_pairs$snp_lead)), .(SNP, donor_id, genotype)]
geno_sub <- as.data.frame(geno_sub)



read_pheno <- function(depot, celltype){
  
  pheno_dir <- "/project/adipose/dynamic_qtl/01.input/phenotype"
  # Read phenotype
  if (depot == "SAT") {
    gpc <- 8
    ppc <- 10
  } else {
    gpc <- 6
    ppc <- 10
  }
  
  file <- paste0(pheno_dir, "/", celltype, "_", depot, "_saige_phenotype.RNA.safe_gene_name.with_BC.gPC", gpc, "_pPC", ppc, ".tsv"
  )
  
  pheno_wide <- fread(file, sep = "\t")
  setDT(pheno_wide)
  return(pheno_wide)
}
all_pheno_genes <- list()

for (ct in celltype) {
  sat_cols <- colnames(read_pheno("SAT", ct))
  vat_cols <- colnames(read_pheno("VAT", ct))
  this_cols <- intersect(sat_cols, vat_cols)
  this_cols <- setdiff(this_cols, "barcode")
  all_pheno_genes[[ct]] <- this_cols
}

common_gene_cols <- Reduce(intersect, all_pheno_genes)
common_gene_cols <- intersect(common_gene_cols, gene_list)

df_plot <- data.frame()
df_stats <- data.frame()

for (i in unique(celltype)){
  phenotype_sat <-read_pheno("SAT",i)
  phenotype_sat_sub <- phenotype_sat %>% select("barcode", all_of(common_gene_cols))
  phenotype_sat_sub$Depot = "SAT"
  phenotype_sat_sub$celltype <- i
  phenotype_vat <-read_pheno("VAT",i)
  phenotype_vat_sub <- phenotype_vat %>% select("barcode", all_of(common_gene_cols))
  phenotype_vat_sub$Depot = "VAT"
  phenotype_vat_sub$celltype <- i
  dat <- rbind(phenotype_sat_sub, phenotype_vat_sub)
  dat <- merge(dat, donor_meta, by = "barcode")
  setDT(dat)
  
  gene_cols <- setdiff(colnames(dat), c("barcode", "Depot", "celltype", "donor_id"))
  gene_cols <- common_gene_cols
  pb <- dat[,
            lapply(.SD, sum, na.rm = TRUE),
            by = .(donor_id, Depot),
            .SDcols = gene_cols
  ]
  pb <- as.data.frame(pb)
  
  pb_meta <- pb[, c("donor_id", "Depot")]
  pb_mat <- as.matrix(pb[, gene_cols, drop = FALSE])
  rownames(pb_mat) <- paste(pb_meta$donor_id, pb_meta$Depot, sep = "__")
  
  dge <- DGEList(counts = t(pb_mat))
  dge <- calcNormFactors(dge)
  
  norm_cpm <- cpm(dge, normalized.lib.sizes = TRUE, log = FALSE)
  log_norm_cpm <- log2(norm_cpm + 1)
  
  pb_norm <- as.data.frame(t(log_norm_cpm))
  pb_norm$donor_id <- pb_meta$donor_id
  pb_norm$Depot <- pb_meta$Depot
  
  pb_geno <- merge(pb_norm, geno_sub, by = "donor_id")
  pb_geno$celltype <- i
  df_plot <- rbind(df_plot, pb_geno)
  # cis stats
  for (g in unique(gene_cols)) {
    
    sat_file <- paste0("/project/01.adipose_sqtl/data/09.eQTL/03.Saige_QTL/02.output/03.formal_calling/ancestry_specific_version/SAT/", i, "/", g, "_cis")
    if (file.exists(sat_file)) {
      lead_sat <- unique(sub$snp_lead[sub$gene_id == g & sub$depot == "SAT"])
      cis_sat <- read.table(sat_file, sep = "\t", header = TRUE)
      cis_sat <- cis_sat[cis_sat$MarkerID %in% lead_sat, , drop = FALSE]
      if (nrow(cis_sat) > 0) {
        cis_sat$Depot <- "SAT"
        cis_sat$celltype <- i
        cis_sat$gene <- g
        df_stats <- rbind(df_stats, cis_sat)
      }
    }
    
    vat_file <- paste0("/project/01.adipose_sqtl/data/09.eQTL/03.Saige_QTL/02.output/03.formal_calling/ancestry_specific_version/VAT/", i, "/", g, "_cis")
    if (file.exists(vat_file)) {
      lead_vat <- unique(sub$snp_lead[sub$gene_id == g & sub$depot == "VAT"])
      cis_vat <- read.table(vat_file, sep = "\t", header = TRUE)
      cis_vat <- cis_vat[cis_vat$MarkerID %in% lead_vat, , drop = FALSE]
      if (nrow(cis_vat) > 0) {
        cis_vat$Depot <- "VAT"
        cis_vat$celltype <- i
        cis_vat$gene <- g
        df_stats <- rbind(df_stats, cis_vat)
      }
    }
  }
  
  
}


## visualize
celltype_order <- c(
  "Early_preadipocyte",
  "Committed_preadipocyte",
  "Adipocyte"
)

library(ggplot2)
library(dplyr)

#------------------------------------------------------------
# 1. 固定 cell type 顺序
#------------------------------------------------------------
celltype_order <- c(
  "Early_preadipocyte",
  "Committed_preadipocyte",
  "Adipocyte"
)

df_plot$celltype <- factor(df_plot$celltype, levels = celltype_order)
#------------------------------------------------------------
# 2. 绘图函数（box + 点，不 jitter）
#------------------------------------------------------------
plot_gene_box <- function(gene, depot, snp, df_plot, outdir) {
  
  df_sub <- df_plot %>%
    filter(
      Depot == depot,
      SNP == snp,
      !is.na(genotype)
    )
  
  if (nrow(df_sub) == 0) return(NULL)
  
  if (!gene %in% colnames(df_sub)) return(NULL)
  
  df_long <- df_sub %>%
    select(donor_id, celltype, genotype, all_of(gene)) %>%
    rename(expr = all_of(gene)) %>%
    filter(!is.na(expr), expr != 0)
  df_long$genotype <- factor(as.character(df_long$genotype), levels = c("0","1","2"))
  if (nrow(df_long) == 0) return(NULL)
  
  p <- ggplot(
  df_long,
  aes(
    x = celltype,
    y = expr,
    fill = factor(genotype),
    group = interaction(celltype, genotype)
  )   ) +
  geom_violin(trim = TRUE, fill = NA, linewidth = 0.8,aes(color = factor(genotype)), width = 0.81,position = position_dodge(width = 0.8)) +
  geom_boxplot(width = 0.16, outlier.shape = NA, fill = "white", color = "black", linewidth = 0.35,
    position = position_dodge(width = 0.8))  +
    theme_bw() +
    labs(
      title = paste0(gene, " (", depot, ", ", snp, ")"),
      x = NULL,
      y = "log2(CPM + 1)",
      fill = "Genotype",
      color = "Genotype"
    ) +
      scale_color_manual(values = c("0"="#4F6D80","1"="#D7CEA3","2"="#B98B86")) +
    scale_fill_manual(values = c("0"="#4F6D80","1"="#D7CEA3","2"="#B98B86")) +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1),
      panel.grid = element_blank()
    )
  
  outfile <- file.path(
    outdir,
    paste0(gene, "_", depot, "_", gsub("[:/]", "_", snp), "_boxplot.pdf")
  )
  
  ggsave(
    filename = outfile,
    plot = p,
    width = 4,
    height = 4
  )
  
  return(p)
}




#------------------------------------------------------------
# 3. 批量画 gene × depot
#------------------------------------------------------------
outdir <- "/project/adipose/dynamic_qtl/02.output/plot"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

for (k in seq_len(nrow(plot_pairs))) {
  
  g <- plot_pairs$gene_id[k]
  d <- plot_pairs$depot[k]
  s <- plot_pairs$snp_lead[k]
  
  message("Plotting: ", g, " ", d, " ", s)
  
  plot_gene_box(
    gene    = g,
    depot   = d,
    snp     = s,
    df_plot = df_plot,
    outdir  = outdir
  )
}







