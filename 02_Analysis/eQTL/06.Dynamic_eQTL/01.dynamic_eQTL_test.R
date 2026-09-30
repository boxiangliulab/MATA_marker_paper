library(dplyr)
library(tidyverse)
library(stringr)
library(data.table)
library(tidyverse)
library(dplyr)
library(glmmTMB)
library(data.table)
library(tidyr)
library(future)
library(future.apply)
library(getopt)
command=matrix(c("depot","d",1,"character",
                 "chr","r",1,"character",
                 "help","h",0,"logical"),byrow=T,ncol=4)
args=getopt(command)

depot=args$depot
chr=args$chr

sat <- read.table("/project/adipose/dynamic_qtl/01.input/SAT_pseudotime_no_eur.tsv",sep = "\t",header = T)
vat <- read.table("/project/adipose/dynamic_qtl/01.input/VAT_pseudotime_no_eur.tsv",sep = "\t",header = T)
meta <- rbind(sat,vat)
####### Prepapre files#######
meta <- meta %>%  # meta: pseudotime_raw, Depot, lineage_annotation
  group_by(Depot) %>%
  mutate(
    pt_raw = slingPseudotime_1,
    mean_early = mean(pt_raw[Annotation_level2 == "Early_preadipocyte"], na.rm = TRUE),
    mean_late  = mean(pt_raw[Annotation_level2 == "Adipocyte"], na.rm = TRUE),
    pt_dir = pt_raw
  ) %>%
  ungroup()
### Scale pt values from 0 to 1
meta <- meta %>%
  group_by(Depot) %>%
  mutate(
    pt01 = (pt_dir - min(pt_dir, na.rm = TRUE)) /
      (max(pt_dir, na.rm = TRUE) - min(pt_dir, na.rm = TRUE))
  ) %>%
  ungroup()
setDT(meta)
meta_small <- meta[, .(barcode, slingPseudotime_1, pt01)]
setkey(meta_small, barcode)
######## Start to test dynamic eQTL #######
library(data.table)

## set file paths
pheno_dir <- "/project/01.adipose_sqtl/data/24.dynamic_qtl/01.dynamic_eQTL/01.input/phenotype/"


# read genotype
geno_raw <- fread("/project/adipose/fullgenotype/Ancestry_specific/00.merge_genotype/genotype_adipose_after_imputation_296_merged_genome_rename_rsID_neur.tsv",sep="\t")
geno_long <- geno_raw %>% pivot_longer(cols = -1, names_to = "donor_id", values_to = "genotype")
geno_long <- data.table(geno_long)

## Read  LD frame
ld_df <- fread("/data/06.genotype/06.plink_file/adipose_after_imputation_296_merged_ld_blocks.ld")

# eQTL files
conditional_cis_eqtl <- read.table("/project/adipose/dynamic_qtl/01.input/SAT_VAT_merged_eQTL.tsv",sep = "\t",header = T)
conditional_cis_eqtl$CHR <- str_split_fixed(conditional_cis_eqtl$top_MarkerID,":",4)[,1]
gene_list <- read.table("/project/adipose/dynamic_qtl/01.input/SAT_VAT_intersect_egene.tsv")
gene_list$V1<- gsub("-","_",gene_list$V1)
conditional_cis_eqtl <- conditional_cis_eqtl[which(conditional_cis_eqtl$gene %in% gene_list$V1),]

prepare_phenotype_for_contex_biased <- function(depot, chr) {
  if (depot == "SAT") {
    pheno_wide <- fread("/project/adipose/dynamic_qtl/01.input/phenotype/SAT_saige_phenotype.RNA.safe_gene_name.with_BC.gPC8_pPC10.tsv", sep = "\t")
  } else {
    pheno_wide <- fread("/project/adipose/dynamic_qtl/01.input/phenotype/VAT_saige_phenotype.RNA.safe_gene_name.with_BC.gPC6_pPC10.tsv", sep = "\t")
  }
  colnames(pheno_wide) <- gsub("-","_",colnames(pheno_wide))
  if (depot == "SAT") {
    meta_cols <- pheno_wide[, 1:26]
    gene_cols <- colnames(pheno_wide)[27:ncol(pheno_wide)]
    gene_cols <- gene_cols[which(gene_cols %in% gene_list$V1)]
  } else {
    meta_cols <- pheno_wide[, 1:24]
    gene_cols <- colnames(pheno_wide)[25:ncol(pheno_wide)]
    gene_cols <- gene_cols[which(gene_cols %in% gene_list$V1)]
  }
  
  cis_qtl_tmp <- conditional_cis_eqtl[which(conditional_cis_eqtl$depot==depot & conditional_cis_eqtl$CHR == chr), ]
  cis_qtl_tmp <- cis_qtl_tmp[which(cis_qtl_tmp$gene %in% gene_cols),]
  cis_qtl_tmp <- as.data.table(cis_qtl_tmp)
  genes_to_keep <- unique(cis_qtl_tmp$gene)
  gene_cols <- intersect(gene_cols, genes_to_keep)
  
  pheno_long <- pheno_wide %>%
    select(all_of(c(colnames(meta_cols), gene_cols))) %>%
    pivot_longer(cols = all_of(gene_cols),
                 names_to = "gene_id",
                 values_to = "expr")
  setDT(pheno_long)
  pheno_long <- meta_small[pheno_long, on = "barcode"]
  setkey(pheno_long, gene_id)
  setDT(geno_long)
  setkey(geno_long, SNP)
  
  lead_snps <- unique(cis_qtl_tmp$top_MarkerID)
  ld_tags <- unique(rbind(
    ld_df[R2 > 0.9 & SNP_A %in% lead_snps, .(lead = SNP_A, tag = SNP_B)],
    ld_df[R2 > 0.9 & SNP_B %in% lead_snps, .(lead = SNP_B, tag = SNP_A)]
  ))
  
  combo_original <- cis_qtl_tmp[, .(gene_id = gene, snp_id = top_MarkerID, snp_lead = top_MarkerID)]
  combo_tag <- merge(cis_qtl_tmp[, .(gene_id = gene, snp_lead = top_MarkerID)],ld_tags, by.x = "snp_lead", by.y = "lead", allow.cartesian = TRUE)[, .(gene_id, snp_id = tag, snp_lead)]
  combo <- unique(rbind(combo_original, combo_tag))
  
  start_time <- Sys.time()
  result <- data.table()
  for (i in 1:nrow(combo)) {
    iter_start <- Sys.time()
    
    gene <- combo$gene_id[i]
    snp  <- combo$snp_id[i]
    
    pheno_sub <- pheno_long[.(gene), on = "gene_id"]
    geno_sub  <- geno_long[.(snp), .(donor_id, genotype)]
    
    if (nrow(pheno_sub) == 0 | nrow(geno_sub) == 0) next
    
    df <- merge(pheno_sub, geno_sub, by = "donor_id")
    df$slingPseudotime_1 <- df$slingPseudotime_1
    
    ### ---------- FIT FULL MODEL ----------
    full <- tryCatch({
      if (depot == "SAT") {
        glmmTMB(
          expr ~ genotype * slingPseudotime_1 + diabetic + sex + BMI + Age + log1p_UMI + percent.mt +
            gPC1 + gPC2 + gPC3 + gPC4 + gPC5 + gPC6 + gPC7 + gPC8 +
            pPC1 + pPC2 + pPC3 + pPC4 + pPC5 + pPC6 + pPC7 + pPC8 + pPC9 + pPC10 +
            (1 | donor_id),
          data = df,
          family = nbinom2(),
          control = glmmTMBControl(parallel = 48)
        )
      } else {
        glmmTMB(
          expr ~ genotype * slingPseudotime_1 + diabetic + sex + BMI + Age + log1p_UMI + percent.mt +
            gPC1 + gPC2 + gPC3 + gPC4 + gPC5 + gPC6 +
            pPC1 + pPC2 + pPC3 + pPC4 + pPC5 + pPC6 + pPC7 + pPC8 + pPC9 + pPC10 +
            (1 | donor_id),
          data = df,
          family = nbinom2(),
          control = glmmTMBControl(parallel = 48)
        )
      }
    }, error = function(e) { return(NULL) })
    
    if (is.null(full)) next
    
    ### ---------- FIT NULL MODEL ----------
    null <- tryCatch({
      if (depot == "SAT") {
        glmmTMB(
          expr ~ slingPseudotime_1 + diabetic + genotype + sex + BMI + Age + log1p_UMI + percent.mt +
            gPC1 + gPC2 + gPC3 + gPC4 + gPC5 + gPC6 + gPC7 + gPC8 +
            pPC1 + pPC2 + pPC3 + pPC4 + pPC5 + pPC6 + pPC7 + pPC8 + pPC9 + pPC10 +
            (1 | donor_id),
          data = df,
          family = nbinom2(),
          control = glmmTMBControl(parallel = 48)
        )
      } else {
        glmmTMB(
          expr ~ slingPseudotime_1 + genotype + sex + BMI + Age + log1p_UMI + percent.mt +
            gPC1 + gPC2 + gPC3 + gPC4 + gPC5 + gPC6 +
            pPC1 + pPC2 + pPC3 + pPC4 + pPC5 + pPC6 + pPC7 + pPC8 + pPC9 + pPC10 +
            (1 | donor_id),
          data = df,
          family = nbinom2(),
          control = glmmTMBControl(parallel = 48)
        )
      }
    }, error = function(e) { return(NULL) })
    
    if (is.null(null)) next
    
    ### ---------- LRT TEST ----------
    p_dynamic <- tryCatch({
      anova(null, full)$"Pr(>Chisq)"[2]
    }, error = function(e) NA)
    
    ### ---------- GET INTERACTION COEF ----------
    coef_full <- coef(summary(full))$cond
    est <- coef_full["genotype:slingPseudotime_1", "Estimate"]
    se  <- coef_full["genotype:slingPseudotime_1", "Std. Error"]
    
    ### ---------- OUTPUT ----------
    tmp <- data.table(
      gene_id = gene,
      snp_id = snp,
      snp_lead = combo$snp_lead[i],
      term = "genotype:slingPseudotime_1",
      p_dynamic = p_dynamic,
      estimate = est,
      se = se,
      depot = depot
    )
    
    result <- rbind(result, tmp)
    
    iter_time <- round(as.numeric(difftime(Sys.time(), iter_start, units = "secs")), 2)
    message(gene, " : ", snp, " done! Time used: ", iter_time, " seconds.")
  }
  #write.table(result,paste0("/data/projects/11003054/e1101943/project/01.adipose_sqtl/data/24.dynamic_qtl/01.dynamic_eQTL/02.output/",depot,"/dynamic_eQTL/",chr,".dynamic_eQTL.tsv"),sep = "\t",quote= F)
  end_time <- Sys.time()
  total_time <- round(difftime(end_time, start_time, units = "mins"), 2)
  message(paste0(depot, " chr", chr, ": ==== Total time used: ", total_time, " minutes ===="))
  return(as.data.frame(result))
}

out <- prepare_phenotype_for_contex_biased(depot,chr)
out_dir <- paste0("/project/adipose/dynamic_qtl/02.output/",depot,"/dynamic_eQTL")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write.table(
  out,
  paste0(out_dir,"/",chr,".dynamic_eQTL.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

