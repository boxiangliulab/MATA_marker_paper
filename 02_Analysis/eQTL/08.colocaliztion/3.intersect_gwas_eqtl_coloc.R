# This script is to generate loci files for colocalization from SAIGE-QTL cis files (We only test eGene)
# Author: Yihan Tong
# Date: 2025.6.24
# Output: cis-eQTL loci with p-value < 1e-3, not in HLA region, not in the region of 500kb of significant SNP

library(data.table)
library(tidyverse)
library(dplyr)
library(glmmTMB)
library(data.table)
library(tidyr)
library(future)
library(future.apply)
library(getopt)
library(stringr)
library("coloc")
command=matrix(c("depot","d",1,"character",
                 "celltype","c",1,"character",
                 "help","h",0,"logical"),byrow=T,ncol=4)
args=getopt(command)

depot=args$depot
celltype=args$celltype

########## Set dir ######
input_dir = paste0("/data/projects/11003054/e1101943/project/01.adipose_sqtl/data/14.coloc/eQTL_loci/",depot,"/1e-3/")
output_inter_dir = paste0("/data/projects/11003054/e1101943/project/01.adipose_sqtl/data/14.coloc/intersect_eQTL_GWAS_loci/",depot,"/")
output_coloc_dir = paste0("/data/projects/11003054/e1101943/project/01.adipose_sqtl/data/14.coloc/colocalization_result/",depot,"/")

######### Read files #####
# MAF of MAATA
MAF <- fread("/data/projects/11003054/e1101943/project/01.adipose_sqtl/data/14.coloc/MAATA_MAF/MAATA_MAF_3pop_296_merged.tsv")
MAF$variant_id <- gsub("^chr", "", MAF$variant_id)
MAF[, MAF_Chinese := Chinese]
MAF[, MAF_Indian  := Indian]
MAF[, MAF_Malay   := Malay]
# one MAF vector for coloc input (minimal-change choice: mean across 3 pops)
MAF[, MAF_mean := rowMeans(.SD, na.rm = TRUE), .SDcols = c("MAF_Chinese","MAF_Indian","MAF_Malay")]

# keep only necessary columns to avoid accidental name collision
MAF <- MAF[, .(variant_id, rsID, variant, MAF_Chinese, MAF_Indian, MAF_Malay, MAF_mean)]

# N of GWAS study
sample_size_gwas <- read.table("/data/projects/11003054/e1101943/project/01.adipose_sqtl/data/14.coloc/GWAS_trait.txt",sep = "\t",header = T)
sample_size_gwas$Study <- gsub("_hg38","",sample_size_gwas$Study)
sample_size_gwas$N <- gsub(",","", sample_size_gwas$N) %>% as.numeric()

# sample_size_gwas <- read.table("/data/projects/11003054/e1101943/project/01.adipose_sqtl/data/10.scDRS/01.magma/01.input/GWAS_trait.txt")
# sample_size_gwas$V1 <- gsub("_hg38","",sample_size_gwas$V1)

filelist <- list.files(
  "/data/projects/11003054/share/data/GWAS_SumStats/prepare_required_columns/threshold_1e-7",
  pattern = "hg38_1e_7\\.txt$",
  full.names = TRUE
)
trait_names <- gsub("_hg38_1e_7\\.txt$", "", basename(filelist))

for (trait in trait_names) {
  message(">>> Processing trait: ", trait)
  
  dir.create(file.path(output_inter_dir, trait), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_coloc_dir, trait), recursive = TRUE, showWarnings = FALSE)
  
  message("  - Loading GWAS summary stats...")
  GWAS_pre <- fread(
    paste0("/data/projects/11003054/share/data/GWAS_SumStats/prepare_required_columns/threshold_1e-7/", trait, "_hg38_1e_7.txt"),
    sep = "\t", header = TRUE
  )
  GWAS_raw <- fread(
   paste0("/data/projects/11003054/share/data/GWAS_SumStats/prepare_required_columns/", trait, "_hg38.txt"),
  
  )  
  if (nrow(GWAS_pre) == 0) {
    message(paste0("  - No trait found in", trait))
    next
  }
  
  GWAS_pre <- GWAS_pre[, .(variant_id = SNPID_hg38, beta = BETA, se = SE, Pvalue = Pvalue)]
  GWAS_raw <- GWAS_raw[, .(variant_id = SNPID_hg38, beta = BETA, se = SE, Pvalue = Pvalue)]
  
  message("  - Loading eQTL nominal loci...")
  nominal <- fread(paste0(input_dir, depot, ".", celltype, ".1e-3.loci"))
  colnames(nominal) <- c("CHR", "POS", "variant_id", "Allele1", "Allele2", "AC_Allele2",
                         "AF_Allele2", "MissingRate", "BETA", "SE", "Tstat", "var", "p.value",
                         "p.value.NA", "Is.SPA", "N", "eGene")
  nominal$variant_id <- gsub("chr", "", nominal$variant_id)
  nominal$variant_id <- sub("(:[^:]+){2}$", "", nominal$variant_id)
  
  common_snp <- merge(nominal, GWAS_pre, by = "variant_id", all = FALSE)
  common_snp <- na.omit(common_snp)
  
  if (nrow(common_snp) == 0) {
    message("  - No intersecting SNPs found. Skipping.")
    next
  }
  
  message("  - Filtering out HLA region...")
  chr_pos <- do.call(rbind, strsplit(common_snp$variant_id, ":"))
  chr <- chr_pos[, 1]
  pos <- as.numeric(chr_pos[, 2])
  hla_filter <- !(chr == "6" & pos >= 25000000 & pos <= 35000000)
  common_snp <- common_snp[hla_filter, ]
  
  if (nrow(common_snp) == 0) {
    message("  - All SNPs filtered by HLA. Skipping.")
    next
  }
  
  message("  - Clumping SNPs...")
  # Perform clumping to retain the most significant SNPs per chromosome,
  # ensuring that no two retained SNPs are within 500 kb of each other.
  # For each chromosome, retain the top SNPs based on p-value,
  # and remove nearby SNPs within ±500 kb to avoid redundant signals.
  common_snp <- as.data.table(common_snp)
  common_snp[, c("chr", "pos") := tstrsplit(variant_id, ":", fixed = TRUE)]
  common_snp[, pos := as.numeric(pos)]
  common_snp[, keep := FALSE]
  
  for (current_chr in unique(common_snp$chr)) {
    dt_chr <- common_snp[chr == current_chr]
    dt_chr <- dt_chr[order(Pvalue)]
    dt_chr[, keep_tmp := TRUE]
    
    for (i in 1:nrow(dt_chr)) {
      if (dt_chr$keep_tmp[i]) {
        current_pos <- dt_chr$pos[i]
        dt_chr[
          pos > current_pos & pos - current_pos < 500000,
          keep_tmp := FALSE
        ]
      }
    }
    common_snp[chr == current_chr & variant_id %in% dt_chr[keep_tmp == TRUE]$variant_id, keep := TRUE]
  }
  
  common_snp_filtered <- common_snp[keep == TRUE]
  message("  - ", nrow(common_snp_filtered), " SNPs retained after clumping.")
  
  write.table(
    common_snp_filtered[, !c("chr", "pos", "keep")],
    paste0(output_inter_dir, trait, "/", depot, ".", celltype, ".ciseQTL_GWAS_inter.loci"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  unique_egene_need_test <- unique(common_snp_filtered$eGene)

  if (nrow(common_snp_filtered) > 0) {
    message("  - Starting coloc.abf for retained SNPs...")
    # test_rs <- data.frame(matrix(NA, nrow = nrow(common_snp_filtered), ncol = 10))
    test_rs <- data.frame(matrix(NA, nrow = length(unique_egene_need_test), ncol = 16))
    #colnames(test_rs) <- c("variant_id", "gene", "nsnps", "PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf","Pvalue_gwas","Pvalue_eqtl")
    colnames(test_rs) <- c("variant_id", "gene", "nsnps","PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf","Pvalue_gwas", "Pvalue_eqtl","beta_gwas", "se_gwas", "beta_eqtl","MAF_Chinese", "MAF_Indian", "MAF_Malay")

    for (j in 1:length(unique_egene_need_test)) {
      egene <- unique_egene_need_test[j]
      coloc_file <- paste0("/data/projects/11003054/e1101943/project/01.adipose_sqtl/data/09.eQTL/03.Saige_QTL/02.output/03.formal_calling/ancestry_specific_version/", depot, "/", celltype, "/", egene, "_cis")
      if (!file.exists(coloc_file)) next
      
      snp_1mb <- fread(coloc_file, sep = "\t", header = TRUE)
      snp_1mb$gene <- egene
      snp_1mb <- snp_1mb[, .(gene, MarkerID, CHR, POS, BETA_eqtl = BETA, p.value_eqtl = p.value, N_eqtl = N)]
      snp_1mb <- snp_1mb[order(p.value_eqtl)]
      snp_1mb[, variant_id := gsub("chr", "", MarkerID)]
      snp_1mb[, variant_id := sapply(strsplit(variant_id, ":"), function(x) paste(x[1], x[2], sep = ":"))]
      top_esnp <- snp_1mb[1,"variant_id"]      
      #test_rs[j, 1:2] <- snp_1mb[1, .(variant_id, gene)]
      #test_rs[j, 10] <- snp_1mb[1, "p.value_eqtl"]
      test_rs[j, c("variant_id","gene")] <- snp_1mb[1, .(variant_id, gene)]
      test_rs[j, "Pvalue_eqtl"] <- snp_1mb[1, p.value_eqtl]
      test_rs[j, "beta_eqtl"]   <- snp_1mb[1, BETA_eqtl]

      input <- merge(snp_1mb, MAF, by = "variant_id")
      GWAS_tmp <- copy(GWAS_raw)
      setnames(GWAS_tmp, c("variant_id", "beta", "se", "Pvalue"), c("variant_id", "beta_gwas", "se_gwas", "Pvalue_gwas"))
      input <- merge(input, GWAS_tmp, by = "variant_id")
      input <- input[order(Pvalue_gwas)]
      input[, varbeta := se_gwas^2]
      unique_input <- input[!duplicated(variant_id)]
      unique_input <- unique_input[!is.na(MAF_Chinese) & MAF_Chinese > 0 & MAF_Chinese < 1]
      # record number of SNPs used for coloc
      test_rs[j, "nsnps"] <- nrow(unique_input)
      # record GWAS stats for the top eSNP
      idx_g <- match(top_esnp, GWAS_raw$variant_id)
      if (!is.na(idx_g)) {
              test_rs[j, "Pvalue_gwas"] <- GWAS_raw$Pvalue[idx_g]
              test_rs[j, "beta_gwas"]   <- GWAS_raw$beta[idx_g]
              test_rs[j, "se_gwas"]     <- GWAS_raw$se[idx_g]
      }
      # record ancestry-stratified MAF for the top eSNP
      idx_m <- match(top_esnp, unique_input$variant_id)
      if (!is.na(idx_m)) {
              test_rs[j, "MAF_Chinese"] <- unique_input$MAF_Chinese[idx_m]
              test_rs[j, "MAF_Indian"]  <- unique_input$MAF_Indian[idx_m]
              test_rs[j, "MAF_Malay"]   <- unique_input$MAF_Malay[idx_m]
              }

      
      if (nrow(unique_input) > 10) {
              if (sample_size_gwas$Case_control[which(sample_size_gwas$Study==trait)] == "Yes") {
                         # If GWAS trait is a case-control study
                         result <- coloc.abf(
                         dataset1 = list(pvalues = unique_input$Pvalue_gwas, type = "cc", beta = unique_input$beta_gwas,
                         varbeta = unique_input$varbeta, s = as.numeric(sample_size_gwas$Case_proportion[sample_size_gwas$Study == trait]),
                         snp = unique_input$variant_id),
                         dataset2 = list(pvalues = unique_input$p.value_eqtl, type = "quant",
                         N = unique(snp_1mb$N_eqtl), snp = unique_input$variant_id),MAF = unique_input$MAF_Chinese
                         )
                      } else if (sample_size_gwas$Case_control[which(sample_size_gwas$Study==trait)] == "No"){
                              # If GWAS trait is not a case-control study
                              result <- coloc.abf(
                              dataset1 = list(pvalues = unique_input$Pvalue_gwas, type = "quant", beta = unique_input$beta_gwas,
                              varbeta = unique_input$varbeta, N = as.numeric(sample_size_gwas$N[sample_size_gwas$Study == trait]),
                              snp = unique_input$variant_id),
                              dataset2 = list(pvalues = unique_input$p.value_eqtl, type = "quant",
                              N = unique(snp_1mb$N_eqtl), snp = unique_input$variant_id),MAF = unique_input$MAF_Chinese
                              )
                              }
        test_rs[j, 3:8] <- t(as.data.frame(result$summary))[1, 1:6]
      }
    }
    
    # extract Pvalue of top variant in eQTL
    test_rs$Pvalue_gwas <- GWAS_raw$Pvalue[match(test_rs$variant_id, GWAS_raw$variant_id)]

    test_rs$depot=depot
    test_rs$celltype=celltype
    test_rs$trait=trait
    write.table(
      test_rs,
      paste0(output_coloc_dir, trait, "/", depot, ".", celltype, "_coloc.tsv"),
      sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE
    )
    message("  - Coloc results written to file.")
  }
  message(">>> Finished trait: ", trait, "\n")
}





