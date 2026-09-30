#### This script is to
### 1. prepare the phenotype for cell-state cis-eQTL
### 2. perform linear regression fitting with poisson distribution
### This analysis will use all conditional cis-eQTL with FDR <= 0.05
### Date: 2025.8.19
### Author: Yihan Tong

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
                 "pc","p",1,"character",
                 "help","h",0,"logical"),byrow=T,ncol=4)
args=getopt(command)

depot=args$depot
chr=args$chr
pc = paste0("pPC", args$pc)

# Project root is 8 levels up from this script
args0 <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args0[grep("--file=", args0)])
project_root <- normalizePath(file.path(dirname(normalizePath(script_path)),
                                        paste(rep("..", 8), collapse = "/")))

outdir <- file.path(project_root, "data/09.eQTL/03.Saige_QTL/02.output/11.cell_state_eQTL/02.adipogenesis", depot, "output")
options(future.globals.maxSize = 2000 * 1024^8)
### Read raw genotype
geno_raw <- fread("/home/project/11003054/e1353486/project/adipose/fullgenotype/Ancestry_specific/00.merge_genotype/adipose_after_imputation_296_merged_genome_rename_rsID_nochr_neur_dosage_header.tsv",sep="\t")
geno_long <- geno_raw %>% pivot_longer(cols = -(1:3), names_to = "donor_id", values_to = "genotype")
geno_long <- data.table(geno_long)

conditional_cis_eqtl <- read.table(file.path(project_root, "data/09.eQTL/03.Saige_QTL/02.output/03.formal_calling/ancestry_specific_version/All_eQTL_SAT_and_VAT_with_beta.tsv"), sep = "\t", header = T)
conditional_cis_eqtl <- conditional_cis_eqtl[which(conditional_cis_eqtl$celltype %in% c("Committed_preadipocyte","Early_preadipocyte") & conditional_cis_eqtl$depot == depot),]
conditional_cis_eqtl$gene <- gsub("-","_",conditional_cis_eqtl$gene)
## Read  LD frame
#ld_df <- fread("/data/projects/11003054/e1101943/project/01.adipose_sqtl/data/06.genotype/ancestry_specific_version/adipose_296_LD.ld")
ld_df <- fread(file.path(project_root, "data/06.genotype/06.plink_file/adipose_after_imputation_296_merged_ld_blocks.ld"))
####################################

prepare_phenotype_for_contex_biased <- function(depot, chr, pc= pc) {
  if (depot == "SAT") {
    pheno_dir <- file.path(project_root, "data/09.eQTL/03.Saige_QTL/02.output/11.cell_state_eQTL/SAT")
    gpc <- 8
    ppc <- 10
  } else {
    pheno_dir <- file.path(project_root, "data/09.eQTL/03.Saige_QTL/02.output/11.cell_state_eQTL/VAT")
    gpc <- 6
    ppc <- 10
  }
  pheno_wide <- fread(file.path(outdir, "celltype_state_phenotype.RNA.safe_name.tsv"), sep = "\t")
  setDT(pheno_wide)
  p <- ncol(pheno_wide)
  N <- if (identical(depot, "SAT")) 27L else 25L
  N <- min(N, p)
  first_names <- names(pheno_wide)[seq_len(N)]
  pc <- intersect(pc, names(pheno_wide))      # keep only columns that actually exist
  meta_names <- unique(c(first_names, pc))
  meta_cols <- pheno_wide[, ..meta_names]
  if (length(pc) == 1L) {
          setnames(meta_cols, old = pc, new = "pPC")
          } else if (length(pc) > 1L) {
                  setnames(meta_cols, old = pc, new = paste0("pPC", seq_along(pc)))
                  }
  #idx_middle <- if (p > (N + 15L)) (N + 1L):(p - 15L) else integer(0)
  if (depot == "SAT") {
          gene_cols <- names(pheno_wide)[27:length(names(pheno_wide))]
          } else if (depot == "VAT") {
                  gene_cols <- names(pheno_wide)[25:length(names(pheno_wide))]
                  }
  keep <- c(first_names, pc, gene_cols)
  pheno_wide <- pheno_wide[, ..keep]
  setnames(pheno_wide, old = pc, new = "pPC")
  cis_qtl_tmp <- conditional_cis_eqtl[which( conditional_cis_eqtl$CHR == chr), ]
  cis_qtl_tmp <- as.data.table(cis_qtl_tmp)
  genes_to_keep <- unique(cis_qtl_tmp$gene)
  gene_cols <- intersect(gene_cols, genes_to_keep)
   pheno_long <- pheno_wide %>%
    select(all_of(c(colnames(meta_cols), gene_cols))) %>%
    pivot_longer(cols = all_of(gene_cols),
                 names_to = "gene_id",
                 values_to = "expr")

  setDT(pheno_long)
  setkey(pheno_long, gene_id)
  setDT(geno_long)
  setkey(geno_long, ID)

  lead_snps <- unique(cis_qtl_tmp$MarkerID)
  ld_tags <- unique(rbind(
    ld_df[R2 > 0.9 & SNP_A %in% lead_snps, .(lead = SNP_A, tag = SNP_B)],
    ld_df[R2 > 0.9 & SNP_B %in% lead_snps, .(lead = SNP_B, tag = SNP_A)]
  ))

  combo_original <- cis_qtl_tmp[, .(gene_id =gene, snp_id = MarkerID, snp_lead = MarkerID)]
  combo_tag <- merge(cis_qtl_tmp[, .(gene_id = gene, snp_lead = MarkerID)],ld_tags, by.x = "snp_lead", by.y = "lead", allow.cartesian = TRUE)[, .(gene_id, snp_id = tag, snp_lead)]
  combo <- unique(rbind(combo_original, combo_tag))

  start_time <- Sys.time()
  result <- data.frame()

  for (i in 1:nrow(combo)) {
    iter_start <- Sys.time()

    gene <- combo$gene_id[i]
    snp <- combo$snp_id[i]

    pheno_sub <- pheno_long[.(gene), ]
    geno_sub <- geno_long[.(snp), .(donor_id, genotype)]

    if (nrow(pheno_sub) == 0 || nrow(geno_sub) == 0) next

    df <- merge(pheno_sub, geno_sub, by = "donor_id")
    df$diabetic <- ifelse(df$Self_defined_DM == "T2DM",1,0)
    df$sex <- ifelse(df$Gender == "F",1,0)
    df$log1p_UMI <- log10(df$nCount_RNA + 1)
    res <- tryCatch({
      if (depot == "SAT") {
        model <- glmmTMB(
          expr ~ genotype * pPC + diabetic + sex + BMI + Age + log1p_UMI + percent.mt +
            gPC1 + gPC2 + gPC3 + gPC4 + gPC5 + gPC6 + gPC7 + gPC8  +
            (1 | donor_id),
          data = df,
           family = nbinom2(),
          control = glmmTMBControl(parallel = 24)
        )
         } else if (depot == "VAT") {
        model <- glmmTMB(
          expr ~ genotype * pPC + diabetic  + sex + BMI + Age + log1p_UMI + percent.mt +
            gPC1 + gPC2 + gPC3 + gPC4 + gPC5 + gPC6 +
            (1 | donor_id),
          data = df,
           family = nbinom2(),
          control = glmmTMBControl(parallel = 24)
        )
      }

      coef <- summary(model)$coefficients$cond
      if (!"genotype:pPC" %in% rownames(coef)) return(NULL)

      row <- coef["genotype:pPC", ]

      # Data frame of interaction term
      tmp_inter <- data.table(
      gene_id = gene,
      snp_id = snp,
      snp_lead = combo$snp_lead[i],
      term = "genotype:pPC",
      estimate = row["Estimate"],
      se = row["Std. Error"],
      pvalue = row["Pr(>|z|)"],
      pPC = pc,
      depot = depot
)
    result <- rbind(tmp_inter, result)
      # Data frame of main effect
      row <- coef["genotype", ]
      tmp_main <- data.table(
      gene_id = gene,
      snp_id = snp,
      snp_lead = combo$snp_lead[i],
      term = "genotype",
      estimate = row["Estimate"],
      se = row["Std. Error"],
      pvalue = row["Pr(>|z|)"],
      pPC = pc,
      depot = depot
)
      result <- rbind(tmp_main, result)

      iter_end <- Sys.time()
      iter_time <- round(as.numeric(difftime(iter_end, iter_start, units = "secs")), 2)
      message(gene, " : ", snp, " done! Time used: ", iter_time, " seconds.")
    }, error = function(e) {
      message("Error in gene ", gene, ", SNP ", snp, ": ", e$message)
      return(NULL)
    })
  }
  end_time <- Sys.time()
  total_time <- round(difftime(end_time, start_time, units = "mins"), 2)
  message(paste0( " chr", chr, "pPC tested",pc,": ==== Total time used: ", total_time, " minutes ===="))
  return(as.data.table(result))
}
    
data <- prepare_phenotype_for_contex_biased(depot=depot,chr=chr,pc = pc)
outdir <- file.path(project_root,
  "data/09.eQTL/03.Saige_QTL/02.output/11.cell_state_eQTL/02.adipogenesis",
  depot, "output/csQTL"
)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# two output filenames
f_intersect <- file.path(outdir, sprintf("%s.%s.chr%s.diabetic_biased.intersect.tsv", depot, pc, chr))
f_main      <- file.path(outdir, sprintf("%s.%s.chr%s.diabetic_biased.main.tsv",      depot, pc, chr))

# if data is a data.table, preferred filtering approach
#fwrite(data[term == "genotype:pPC"], f_intersect, sep = "\t", quote = FALSE)
#fwrite(data[term == "genotype"],     f_main,      sep = "\t", quote = FALSE)

data <- as.data.frame(data)
fwrite(subset(data, term == "genotype:pPC"), f_intersect, sep="\t", quote=FALSE)
fwrite(subset(data, term == "genotype"),     f_main,      sep="\t", quote=FALSE)


