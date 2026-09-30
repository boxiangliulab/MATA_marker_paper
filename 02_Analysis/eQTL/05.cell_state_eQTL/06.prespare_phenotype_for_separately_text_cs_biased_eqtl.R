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
                 "pc","p",1,"character",
                 "help","h",0,"logical"),byrow=T,ncol=4)
args=getopt(command)

depot=args$depot
pc = args$pc
pc_id = pc
pc = paste0("pPC",pc)

# Project root is 8 levels up from this script
args0 <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args0[grep("--file=", args0)])
project_root <- normalizePath(file.path(dirname(normalizePath(script_path)),
                                        paste(rep("..", 8), collapse = "/")))

## significant lead eQTL
sig_context_biased_df <- read.table(file.path(project_root, "data/09.eQTL/03.Saige_QTL/02.output/03.formal_calling/ancestry_specific_version/All_eQTL_SAT_and_VAT_with_beta.tsv"), sep = "\t", header = T)
sig_context_biased_df <- sig_context_biased_df[which(sig_context_biased_df$celltype %in% c("Committed_preadipocyte","Early_preadipocyte") & sig_context_biased_df$depot == depot),]
sig_context_biased_df$gene <- gsub("-","_",sig_context_biased_df$gene)


output <- file.path(project_root, "data/09.eQTL/03.Saige_QTL/02.output/11.cell_state_eQTL/02.adipogenesis", depot, "01.phenotype_for_separate_run")
dir.create(output, recursive = TRUE, showWarnings = FALSE)

# Read phenotype
pheno_dir <- file.path(project_root, "data/09.eQTL/03.Saige_QTL/02.output/11.cell_state_eQTL/02.adipogenesis", depot, "output/celltype_state_phenotype.RNA.safe_name.tsv")
pheno_wide <- fread(pheno_dir, sep = "\t")


all_cols <- colnames(pheno_wide)

keep_cols <- all_cols[
  !grepl("^pPC", all_cols) | all_cols == pc
]

dt <- copy(pheno_wide[, ..keep_cols])
setnames(dt, old = pc, new = "pPC")


if (depot == "SAT") {
        gene_cols <- colnames(dt)[21:ncol(dt)]
        } else if (depot == "VAT"){
                gene_cols <- colnames(dt)[19:ncol(dt)]
                }
  

## =========================
## add log10(UMI + 1)
## =========================
dt[, log10_UMI := log10(nCount_RNA + 1)]
## =========================
## recode Gender
## M = 0, F = 1
## =========================
dt[, Gender := fifelse(Gender == "F", 1L,
                fifelse(Gender == "M", 0L, NA_integer_))]

## =========================
## recode Self_defined_DM
## T2DM = 1, others = 0
## =========================
dt[, Self_defined_DM := fifelse(Self_defined_DM == "T2DM", 1L, 0L)]

## =========================
## make pPC_quantile: 0 / 1 / 2
## by tertile
## =========================
qs <- quantile(dt$pPC, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
dt[, pPC_quantile := cut(
  pPC,
  breaks = qs,
  labels = c(0L, 1L, 2L),
  include.lowest = TRUE
)]

dt[, pPC_quantile := as.integer(as.character(pPC_quantile))]

dt[, barcode := NULL]
dt[, nCount_RNA := NULL]
dt[, pPC := NULL]
dt[, Annotation_level2 := NULL]

meta_cols <- setdiff(names(dt),gene_cols)
dt <- dt[, c(meta_cols, gene_cols), with = FALSE]

## =========================
## split and save
## 0 = low
## 1 = median
## 2 = high
## =========================
suffix_map <- c("low", "median", "high")

for (i in 0:2) {
  sub_dt <- dt[pPC_quantile == i]
  
  outfile <- file.path(
    output,
    paste0("phenotype_", suffix_map[i + 1], ".tsv")
  )
  
  fwrite(sub_dt, outfile, sep = "\t", quote = FALSE, na = "NA")
}







