#############################
##### Author:Yihan Tong #####
##### Date: 2024.12.17 ######
######## If you are using output from Cellranger ##############

.libPaths("./.conda/envs/seurat/lib/R/library")
library(getopt)
command=matrix(c("sample","s",1,"character",
                 "batch","b",1,"character",
                 "help","h",0,"logical"),byrow=T,ncol=4)
args=getopt(command)
if( !is.null(args$help) || is.null(args$sample) )
{ cat(paste(getopt(command,usage=T),"\n"))
  q() }

sample=args$sample
batch=args$batch

######################################################
################## Preparation #######################
######################################################

# Libraries
library(DropletUtils)
library(decontX)
library(gridExtra)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(DoubletFinder)
library(Seurat)
library(tidyverse)
library(stringr)
library(sctransform)
library(Azimuth)
library(randomcoloR)
library(SingleCellExperiment)
library(Azimuth)
library(umap)
library(DoubletFinder)
# Set colors
myColors<-distinctColorPalette(20)
colors <- c("Adipocyte" = "#E35663",
            "ASPC" = "#A590A7",
            "B" = "#E4C551",
            "Dendritic" = "#CAE7E6",
            "Endothelial" = "#B7C3E3",
            "Lymphatic Endothelial" = "#EAC6E1",
            "Macrophage" = "#BAE05E",
            "Mesothelial" = "#629687",
            "Monocyte" = "#E5D198",
            "NK" = "#854891",
            "Pericyte" = "#646BDE",
            "Smooth Muscle" = "#E791BF",
            "T" = "#DD89E4",
            "Mast" = "#8537E3",
            "Endometrial" = "#E49151",
            "Neutrophil" = "#8B9447",
            "SAT" = "#33a02c",
            "VAT" = "#ff7f00",
            "Male" = "#1f78b4",
            "Female" = "#e31a1c")

#################### Step 01: Identify empty cells ####################
setwd("./06.seurat_qc/")
if(!dir.exists(batch)){dir.create(batch)}
setwd(batch)
if(!dir.exists(sample)){dir.create(sample)}
setwd(sample)
if(!dir.exists("step01")){dir.create("step01")}
setwd("step01")

# This step is for setting background for removing ambient RNA
set.seed(123)
# First choice: Self-detemine the empty droplets with nUMIs < 500
message("Starting to identify empty cells......")
counts <- Read10X(paste0("./04.cellranger/",batch,"/",sample,"/outs/raw_feature_bc_matrix/"))
seuratObject_raw <- CreateSeuratObject(counts)
seuratObject_empty <- subset(seuratObject_raw,nCount_RNA<=500) # Assume cells with nUMI <= 500 as empty cells
counts.empty <- GetAssayData(object = seuratObject_empty, layer = "counts")
sce_empty <- SingleCellExperiment(list(counts = counts.empty))
print("Step01 done!Congrat!😊")

# Second choice: You can also use emptyDrops from DropletUtils
#message("01. Starting to identify empty cells......")
#my.counts <- Read10X(paste0("./",batch,"/04.cellranger/",sample,"/outs/raw_feature_bc_matrix/"))
#e.out <- emptyDrops(my.counts)
#is.cell <- e.out$FDR <= 0.01
#sce_cells <- counts[, which(is.cell)]
#sce_empty <- counts[, which(!is.cell)]
#message(paste0(c(' ...found ', sum(is.cell), ' cells and ', sum(!is.cell), ' empties')))
#print("Step01 done!Congrat!😊")


#################### Step 02: Ambient RNA correction ####################
setwd(paste0("./06.seurat_qc/",batch,"/",sample))
set.seed(123)
#if(!dir.exists(sample)){dir.create(sample)}
#setwd(sample)
if(!dir.exists("step02")){dir.create("step02")}
setwd("step02")

message("02. Starting to identify ambient RNA......")
# If you want to start with CellRanger filtered count matrix:
counts_filtered <- Read10X(paste0("./04.cellranger/",batch,"/",sample,"/outs/filtered_feature_bc_matrix"))
message(paste0("There are ",ncol(counts_filtered)," nuclei in total identify by CellRanger"))
# Create a SingleCellExperiment object and run decontX
X_sce <- SingleCellExperiment(list(counts = counts_filtered)) # use filtered matrix from cellranger and use empty cells identified last step
set.seed(123)
sce <- decontX(X_sce,background = sce_empty,estimateDelta = TRUE, dbscanEps = 0.5,delta = c(10,5))

#If you want to identify empty cells yourselves:
#X_sce <- SingleCellExperiment(list(counts = sce_cells))
#X_empty <- SingleCellExperiment(list(counts = sce_empty))
#sce <- decontX(X_sce,background = X_empty)


# Create a Seurat object from a SCE with decontX results
options(Seurat.object.assay.calcn = TRUE)
adipose_Seurat <- CreateSeuratObject(round(decontXcounts(sce)),min.cells = 3,min.features  = 10,project = "adipose",assay = "RNA")
#adipose_Seurat[["percent.mt"]] <- PercentageFeatureSet(adipose_Seurat, pattern = "^MT-")
adipose_Seurat <- PercentageFeatureSet(adipose_Seurat, pattern = "^MT-", col.name = "percent.mt")
adipose_Seurat <- PercentageFeatureSet(adipose_Seurat, pattern = "^[Rr][Pp][LlSs]", col.name = "percent.ribo")
# Add contamination ratio in the metadata
cont_ratio <- as.data.frame(sce@metadata$decontX$estimates$all_cells$contamination)
names(cont_ratio) <- "contamination"
adipose_Seurat <- AddMetaData(adipose_Seurat,cont_ratio,"contamination")
print("Step02 done!Congrat!😊")
adipose_Seurat@meta.data$barcode <- rownames(adipose_Seurat@meta.data)
write.table(adipose_Seurat@meta.data[,c("barcode","nCount_RNA","nFeature_RNA","percent.mt")],"02.metadata_no_filter_demuxlet.tsv",row.names = F, quote = F, sep = "\t")
#################### Step 03: Add individual in metadata ####################
demuxlet <- read.table(paste0("./02.demultiplex/",batch,"/",sample,"/",sample,".best"),sep = "\t",header = T)
original_cell_num <- nrow(demuxlet)
demuxlet$BARCODE <- gsub("\\-1",paste0("_",sample),demuxlet$BARCODE)
demuxlet$cell <-  demuxlet$BARCODE
demuxlet$batch <- sample
demuxlet$donor_id <- "a"
demuxlet$donor_id <- ifelse(demuxlet$DROPLET.TYPE=="DBL","Doublet",demuxlet$SNG.BEST.GUESS)
table(demuxlet$donor_id)
demultiplex <- demuxlet
demultiplex <- column_to_rownames(demultiplex,var = "BARCODE")
demultiplex <- demultiplex[rownames(adipose_Seurat@meta.data),]
adipose_Seurat <- AddMetaData(adipose_Seurat,metadata = demultiplex)
rownames(adipose_Seurat@meta.data) <- adipose_Seurat@meta.data$barcode

#################### Step 04: Visualize initial QC matrixes ####################
setwd(paste0("./06.seurat_qc/",batch,"/",sample))
set.seed(123)
#if(!dir.exists(sample)){dir.create(sample)}
#setwd(sample)
if(!dir.exists("step03")){dir.create("step03")}
setwd("step03")
if(!dir.exists("demuxlet")){dir.create("demuxlet")}
setwd("demuxlet")
message("03. Starting to visualize......")
all_ngene <- data.frame(nGene = adipose_Seurat$nFeature_RNA,
                        nUMI = adipose_Seurat$nCount_RNA,
                        pMito = adipose_Seurat$percent.mt)
nGene_tmp <- as.data.frame(adipose_Seurat@meta.data)

# Assign Cell-Cycle Scores
Idents(adipose_Seurat) <- "barcode"
adipose_Seurat <- Seurat::NormalizeData(adipose_Seurat)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
adipose_Seurat <- CellCycleScoring(adipose_Seurat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
DefaultAssay(adipose_Seurat) <- 'RNA'
## Reduction + cluster
options(future.globals.maxSize = 8000 * 1024^2)
adipose_Seurat <- SCTransform(adipose_Seurat,verbose = FALSE,vars.to.regress = c("percent.mt", "percent.ribo","G2M.Score","S.Score"),vst.flavor = "v2")
adipose_Seurat <- Seurat::RunPCA(adipose_Seurat, assay="SCT", verbose = FALSE)
adipose_Seurat <- Seurat::RunUMAP(adipose_Seurat, reduction = "pca",
                       dims = 1:30, verbose = FALSE)
adipose_Seurat <- Seurat::FindNeighbors(adipose_Seurat, reduction = "pca",
                             dims = 1:30, verbose = FALSE)

adipose_Seurat <- Seurat::FindClusters(
    object = adipose_Seurat
)
print("Finish SCT!")
# Run automatic annotation
DefaultAssay(adipose_Seurat) <- 'RNA'
adipose_Seurat <- RunAzimuth(adipose_Seurat, reference = "./01.Azimth_adipose_ref/adiposeref.SeuratData/inst/azimuth/")
print("01.finish Reading raw data") # change the reference, or find a way to annotate yourself, this is for doublet identification

p1 <- DimPlot(adipose_Seurat, reduction = "ref.umap",group.by = "predicted.celltype.l1",cols = myColors[1:length(unique(adipose_Seurat$predicted.celltype.l1))]) + ggtitle("Raw cells before any QC")
p2 <- DimPlot(adipose_Seurat, reduction = "umap",group.by = "predicted.celltype.l1",cols = myColors[1:length(unique(adipose_Seurat$predicted.celltype.l1))]) + ggtitle("Raw cells before any QC")
pdf("02.annotation_result_before_removing_doublet_SCT.pdf",width = 8,height = 10)
p1/p2
dev.off()

## About metadata distribution
# Individual distribution
individual <- as.data.frame(table(adipose_Seurat@meta.data$donor_id))
individual_col <- c()
for (i in 1:16) {
    ind <- as.character(individual$Var1[!individual$Var1 %in% c("unassigned", "doublet")])[i]
    cols <- c("#3cb44b","#ffe119","#4363d8","#f58231","#911eb4","#46f0f0","#f032e6","#bcf60c","#fabebe","#008080","#e6beff","#9a6324","#fffac8","#800000","#aaffc3","#808000")[i]
    names(cols) <- ind
    individual_col <- c(cols,individual_col)
}
individual_col <- c(individual_col,"unassigned" = "#bababa","doublet" = "#1a1a1a")



write.table(individual,"individual_distribution_before_doublet_removel.tsv",row.names=F,quote=F,sep="\t")
print("02.finish plot individual")

singlets_num <- sum(individual$Freq)-individual$Freq[which(individual$Var1 == "doublet")]
individual$Var1 <- factor(individual$Var1,levels=individual$Var1[order(individual$Freq)])
y_intercept <- round(singlets_num / 16)
p1 <- ggplot(individual, aes(x = Var1, y = Freq, fill = Var1)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = individual_col) +
    theme_minimal() +
    labs(title = "Bar Plot of Frequencies", x = "Individual", y = "Count") +
    geom_hline(yintercept = y_intercept, linetype = "dashed", color = "gray") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

p2 <- DimPlot(adipose_Seurat, reduction = "umap",group.by = "donor_id",cols=individual_col)
pdf("03.individual_distribution_before_debris_removel.pdf",width=5,height=8)
p1/p2
dev.off()


# Plot XIST gene expression
new_strings <- as.character(individual$Var1[!individual$Var1 %in% c("unassigned", "doublet")])
adipose_Seurat_filter<-subset(adipose_Seurat,cells=rownames(adipose_Seurat@meta.data)[which(adipose_Seurat$donor_id %in% new_strings)])
pdf("04.sex_identity_before_removing_doublet.pdf",width=40,height=4)
FeaturePlot(adipose_Seurat_filter, features = "XIST",split.by = "donor_id")
dev.off()

# Plot celltype contamination dist√ribution
cont_all <- ggplot(as.data.frame(adipose_Seurat@meta.data),aes(x=contamination))+
        geom_histogram(binwidth = 0.1, fill = "grey", color = "black")
        ggtitle("Contamination rate of nuclei before doublet removel")+
        theme_classic()
cont_ct <- ggplot(as.data.frame(adipose_Seurat@meta.data),aes(x=contamination))+
        geom_histogram(binwidth = 0.05, fill = "grey", color = "black")+
        facet_wrap(~ predicted.celltype.l1, ncol = 3)+
        ggtitle("Contamination rate of each celltype before doublet removel")+
        theme_classic()
pdf("05.contamination_rate_b4_DBL_removal.pdf",width=5,height=10)
cont_all/cont_ct
dev.off()
print("Step03 done!Congrat!😊")

adipose_Seurat <- subset(adipose_Seurat,subset = nFeature_RNA >300 &
                           nCount_RNA > 500 &
                           nCount_RNA < 200000  &
                            percent.mt <= 5)

message(paste0("There are ", nrow(adipose_Seurat@meta.data)," nuclei left after remove according to the nFeature, mt and count"))
p1 <- DimPlot(adipose_Seurat, reduction = "ref.umap",group.by = "predicted.celltype.l1",cols = myColors[1:length(unique(adipose_Seurat$predicted.celltype.l1))]) + ggtitle("Cells after remove low count")
p2 <- DimPlot(adipose_Seurat, reduction = "umap",group.by = "predicted.celltype.l1",cols = myColors[1:length(unique(adipose_Seurat$predicted.celltype.l1))]) + ggtitle("Cells after remove low count")
pdf("06.annotation_result_after_removing_low_count_SCT.pdf",width = 8,height = 10)
p1/p2
dev.off()

#message(paste0("There are ", nrow(adipose_Seurat@meta.data)," nuclei left after remove unassigned"))
DefaultAssay(adipose_Seurat) <- 'SCT'

setwd(paste0("./06.seurat_qc/",batch,"/",sample))
set.seed(123)
#if(!dir.exists(sample)){dir.create(sample)}
#setwd(sample)
if(!dir.exists("seurat_rds")){dir.create("seurat_rds")}
setwd("seurat_rds")
saveRDS(adipose_Seurat,"seurat_object_after_filter_low_gene_before_DBL_removel_SCT_demuxlet.rds")
# Stop here, you can check the metadata QC plot. Indentify the sex & fat type youself and add this info in the metadata.





