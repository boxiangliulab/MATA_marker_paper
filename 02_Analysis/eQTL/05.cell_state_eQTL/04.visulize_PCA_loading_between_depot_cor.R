library(tidyverse)
library(linkET)
library(Hmisc)
library(corrplot)

# Project root is 8 levels up from this script
args0 <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args0[grep("--file=", args0)])
project_root <- normalizePath(file.path(dirname(normalizePath(script_path)),
                                        paste(rep("..", 8), collapse = "/")))

base_data <- file.path(project_root, "data/09.eQTL/03.Saige_QTL/02.output/11.cell_state_eQTL/02.adipogenesis")

sat <- read.table(file.path(base_data, "SAT/output/SAT_adipogenesis_PCA_loadings_PC1-8.tsv"), sep="\t", header=T)
vat <- read.table(file.path(base_data, "VAT/output/VAT_adipogenesis_PCA_loadings_PC1-8.tsv"), sep="\t", header=T)

common_gene <- intersect(sat$gene, vat$gene)

sat <- sat[sat$gene %in% common_gene, ]
vat <- vat[vat$gene %in% common_gene, ]

sat <- sat[match(common_gene, sat$gene), ]
vat <- vat[match(common_gene, vat$gene), ]

############################
## 4. Extract PC matrix
############################

sat_mat <- as.matrix(sat[, -1])
vat_mat <- as.matrix(vat[, -1])

############################
## 5. Compute Spearman correlation
## VAT PCs vs SAT PCs
############################

cor_cross <- cor(
  vat_mat,
  sat_mat,
  method = "spearman",
  use = "pairwise.complete.obs"
)

############################
## 6. Rename axes
############################

colnames(cor_cross) <- paste0("SAT_PC_",1:ncol(cor_cross))
rownames(cor_cross) <- paste0("VAT_PC_",1:nrow(cor_cross))

############################
## 7. Define color palette
############################

pal11 <- c(
  "#8e0152","#c51b7d","#de77ae","#f1b6da","#fde0ef",
  "#f7f7f7",
  "#e6f5d0","#b8e186","#7fbc41","#4d9221","#276419"
)
############################
## 8. Plot correlation
############################

write.table(
  as.data.frame(cor_cross),
  file.path(base_data, "SAT_VAT_PCA_loading_cor.tsv"),
  sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA
)

pdf(file.path(base_data, "SAT_VAT_PCA_loading_cor.pdf"), width=5, height=5)
corrplot(
  cor_cross,
  type = "lower",
  method = "square",
  col = pal11,
  cl.lim = c(-1,1),
  border = "black",
  diag = FALSE,
  tl.col = "black",
  tl.srt = 45,
  tl.cex = 1.2
)
dev.off()
