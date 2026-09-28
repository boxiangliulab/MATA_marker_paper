# last modified: 2025.02.11
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
#################### Step 04: Doublet Identify for each library ####################
# Here we emplyed DoubletFinder, you can also use others
# Try to find doublet by each individual, but it returns the exactly the same result as treating a library as a whole
setwd(paste0("./06.seurat_qc/",batch,"/",sample))
set.seed(123)
#if(!dir.exists(sample)){dir.create(sample)}
#setwd(sample)
if(!dir.exists("step04")){dir.create("step04")}
setwd("step04")
if(!dir.exists("demuxlet")){dir.create("demuxlet")}
setwd("demuxlet")
message("04. Starting to identify doublets......")


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
myColors<-distinctColorPalette(40)
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

# Add library info in the metadata
adipose_Seurat <- readRDS(paste0("./06.seurat_qc/",batch,"/",sample,"/seurat_rds/seurat_object_after_filter_low_gene_before_DBL_removel_SCT_demuxlet.rds"))
adipose_Seurat$batch <- sample
adipose_Seurat$group <- batch

# Add metadata (If you haven't done so in first stepQC)
metadata <- as.data.frame(adipose_Seurat@meta.data)
grouplist <- read.table(paste0("./06.seurat_qc/00.meta_manually_identify/Grouplist/",sample,".tsv"),sep = "\t",header=T)
metadata <- merge(metadata,grouplist,by.x = "donor_id",by.y = "TR.ID")
rownames(metadata) <- metadata$barcode
metadata <- metadata[,setdiff(names(metadata),names(adipose_Seurat@meta.data))]
adipose_Seurat <- AddMetaData(adipose_Seurat,metadata)

# Set color for individual
individual <- as.data.frame(table(adipose_Seurat@meta.data$donor_id))
individual_col <- c()
for (i in 1:16) {
    ind <- as.character(individual$Var1[!individual$Var1 %in% c("unassigned", "doublet")])[i]
    cols <- c("#3cb44b","#ffe119","#4363d8","#f58231","#911eb4","#46f0f0","#f032e6","#bcf60c","#fabebe","#008080","#e6beff","#9a6324","#fffac8","#800000","#aaffc3","#808000")[i]
    names(cols) <- ind
    individual_col <- c(cols,individual_col)
}
individual_col <- c(individual_col,"unassigned" = "#bababa","doublet" = "#1a1a1a")

# Read demuxlet result
demuxlet <- read.table(paste0("./02.demultiplex/",batch,"/",sample,"/",sample,".best"),sep = "\t",header = T)
original_cell_num <- nrow(demuxlet)
demuxlet$cell <- gsub("\\-1",paste0("_",sample),demuxlet$BARCODE)
demuxlet$batch <- sample
table(demuxlet$DROPLET.TYPE)
demuxlet$donor_id <- "a"
demuxlet$donor_id <- ifelse(demuxlet$DROPLET.TYPE %in% c("DBL","AMB"),"Doublet",demuxlet$SNG.BEST.GUESS)
demultiplex <-  demuxlet
demultiplex$DROPLET.TYPE <- ifelse(demultiplex$DROPLET.TYPE %in% c("DBL","AMB"), "Doublet", "Singlet")
adipose_Seurat$DROPLET.TYPE <- ifelse(adipose_Seurat$DROPLET.TYPE %in% c("DBL","AMB"), "Doublet", "Singlet")
genetic_DBL_num <- length(which(demultiplex$DROPLET.TYPE == "Doublet"))
DF_doublet <- data.frame()
for (i in unique(adipose_Seurat$donor_id)) {
  if (length(which(adipose_Seurat$donor_id == i)) >= 50) { # we remove those donors with less 50 nuclei in this library
    test <- subset(adipose_Seurat,donor_id==i)
    DefaultAssay(test) <- 'RNA'
    options(future.globals.maxSize = 8000 * 1024^2)
    test <- SCTransform(test,verbose = FALSE,vars.to.regress = c("percent.mt", "percent.ribo","G2M.Score","S.Score"),vst.flavor = "v2")
    test <- Seurat::RunPCA(test, assay="SCT", verbose = FALSE)
    test <- Seurat::RunUMAP(test, reduction = "pca",
                            dims = 1:30, verbose = FALSE)
    test <- Seurat::FindNeighbors(test, reduction = "pca",
                                  dims = 1:30, verbose = FALSE)
    test <- Seurat::FindClusters(
      object = test)
    estimated.doublet.ratio <- (nrow(test@meta.data)/10000)*0.04 # DBL rate if you do not use ground truth
    test <- subset(test,cells=rownames(test@meta.data)[which(test$DROPLET.TYPE %in% c("Doublet","Singlet"))])
    ## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
    sweep.res.list_adipose <- paramSweep(test, PCs = 1:30, sct = TRUE)
    sweep.stats_adipose <- summarizeSweep(sweep.res.list_adipose, GT = FALSE)
    bcmvn_adipose <- find.pK(sweep.stats_adipose)
    pdf(paste0("04.bcmvn_",i,".pdf"),width = 12,height = 5)
    p<-ggplot(bcmvn_adipose,aes(pK,BCmetric,group=1))+geom_point()+geom_line()
    print(p)
    dev.off()
    ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
    homotypic.prop <- modelHomotypic(test$predicted.celltype.l1)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
    nExp_poi <- round(estimated.doublet.ratio*nrow(test@meta.data))
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
    test <- doubletFinder(test, PCs = 1:30, pN = 0.25, pK = as.numeric(as.character(bcmvn_adipose$pK[which.max(bcmvn_adipose$BCmetric)])), nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
    colnames(test@meta.data)[grep("DF.classifications",colnames(test@meta.data))] <- "DF.classfication1"
    test <- doubletFinder(test, PCs = 1:30, pN = 0.25, pK = as.numeric(as.character(bcmvn_adipose$pK[which.max(bcmvn_adipose$BCmetric)])), nExp = nExp_poi.adj, reuse.pANN = "DF.classfication1", sct = TRUE)
    colnames(test@meta.data)[grep("DF.classifications_",colnames(test@meta.data))] <- "DF.classfication2"
    print(paste0(i," has ", length(which(test$DF.classfication2=="Doublet"))," DF doublets"))
    colnames(test@meta.data)[grep("pANN",colnames(test@meta.data))] <- "pANN_value"
    metadata <- test@meta.data[,c("barcode","DF.classfication2","pANN_value")]
    DF_doublet <- rbind(metadata,DF_doublet)

  }
}
print(paste0("There are ", nrow(DF_doublet), " nuclei tested in doubletfinder"))
rownames(DF_doublet) <- DF_doublet$barcode
adipose_Seurat <- AddMetaData(adipose_Seurat,DF_doublet)

# Print doublets found by DF in individual
test <- adipose_Seurat@meta.data[which(adipose_Seurat@meta.data$DF.classfication2=="Doublet"),]
data <- data.frame(table(test$donor_id))
write.table(data,"doublet_found_by_DF_in_each_donor.tsv",sep = "\t",quote=F)

# Plot doublet finder result
doublet.df <- FetchData(adipose_Seurat, vars = c("PC_1","PC_2","DF.classfication2"))
doublet.df <- data.frame(adipose_Seurat@reductions$umap@cell.embeddings) %>% mutate(DF =  adipose_Seurat$DF.classfication2)
dbl_rate<- as.numeric(table(doublet.df$DF)[1]/original_cell_num)

p4 <- ggplot(data = doublet.df,aes(x=umap_1,y =umap_2,color = DF))+geom_point(size = 0.8) +
  theme_classic() + scale_color_manual(values = c("Singlet" = "#1f78b4", "Doublet" = "#e31a1c"))+
  ggtitle(paste0("n = ", nrow(adipose_Seurat@meta.data), ", double rate (DF) = ", round(dbl_rate,digits = 3))) +
  guides(color=guide_legend(override.aes=list(size=4)))

# final singlet after removing doublet finder
final.singlet.barcode <- rownames(doublet.df[which(doublet.df$DF=="Singlet"),])
final.singlet.barcode <- final.singlet.barcode[which(final.singlet.barcode %in% demultiplex$cell[which(demultiplex$DROPLET.TYPE == "Singlet")])] # intersect with result from demultiplex
doublet_final <-  data.frame(adipose_Seurat@reductions$umap@cell.embeddings)
doublet_final$barcode <- rownames(doublet_final)
doublet_final$final.doublet <- "Doublet"
doublet_final[which(doublet_final$barcode %in% final.singlet.barcode),"final.doublet"] <- "Singlet"
print(table(doublet_final$final.doublet))
write.table(doublet_final,"final_doublet.tsv",sep = "\t",quote=F)

# plot result
dbl_rate_2<- as.numeric(table(doublet_final$final.doublet)[1]/original_cell_num)

p5 <- ggplot(data = doublet_final,aes(x=umap_1,y =umap_2,color = final.doublet ))+geom_point(size = 0.8) +
  theme_classic() + scale_color_manual(values = c("Singlet" = "#1f78b4", "Doublet" = "#e31a1c"))+
  ggtitle(paste0("n = ", nrow(adipose_Seurat@meta.data), ", double rate (ALL) = ", round(dbl_rate_2,digits = 3)))+
  guides(color=guide_legend(override.aes=list(size=4)))
#dev.off()
pdf("05.UMAP_of_seurat_for_final_singlets_and_doublets.pdf",width = 12,height = 5)
p4 | p5
dev.off()

colnames(adipose_Seurat@meta.data)[grep("pANN",colnames(adipose_Seurat@meta.data))] <- "pANN_value"

pdf("06.density_plot_of_pANN_between_all_and_singlet_only.pdf",width=4,height=4)
pANN_filtered <- na.omit(adipose_Seurat$pANN_value)
singlet_filtered <- na.omit(adipose_Seurat$pANN_value[final.singlet.barcode])
plot(density(pANN_filtered), ylim = c(0, 30), main = "")
lines(density(singlet_filtered), col = "red")
dev.off()

## hightlighting doublets identified by different tools
doublet.df <- FetchData(adipose_Seurat, vars = c("umap_1","umap_2","DF.classfication2"))
doublet.df$barcode <- rownames(doublet.df)
doublet.df$demultiplex <- "Doublet"
doublet.df[doublet.df$barcode %in% demultiplex$cell[which(demultiplex$DROPLET.TYP == "Singlet")],"demultiplex"] <- "Singlet"
doublet.df$tools_indentify <- ifelse(doublet.df$DF.classfication2 == "Singlet" & doublet.df$demultiplex == "Doublet",
                                    "DF_SNG + demultiplex_DBL",
                                    ifelse(doublet.df$DF.classfication2 == "Singlet" & doublet.df$demultiplex == "Singlet",
                                           "DF_SNG + demultiplex_SNG",
                                           ifelse(doublet.df$DF.classfication2 == "Doublet" & doublet.df$demultiplex == "Doublet",
                                                  "DF_DBL + demultiplex_DBL",
                                                  "DF_DBL + demultiplex_SNG")))


pdf("07.UMAP_of_seurat_for_highlighting_doublets_from_different_tools.pdf",width = 6,height = 4)
ggplot(data = doublet.df,aes(x=umap_1,y =umap_2,color = tools_indentify ))+geom_point(size = 0.3) +
  theme_classic() + scale_color_manual(values = c("DF_SNG + demultiplex_SNG" = "grey", "DF_DBL + demultiplex_SNG" = "red",
                                                  "DF_SNG + demultiplex_DBL" = "blue", "DF_DBL + demultiplex_DBL" = "green"))+
  ggtitle(paste0("n = ", nrow(adipose_Seurat@meta.data), ", double rate (ALL) = ", round(dbl_rate_2,digits = 3))) +
  guides(color=guide_legend(override.aes=list(size=4)))
dev.off()

######## subset final high quality cell ########
high_quality <- c(final.singlet.barcode)
adipose_Seurat <- subset(adipose_Seurat,cells = high_quality)
print(paste0("After remove DBL, there are ", nrow(adipose_Seurat@meta.data)," nuclei left"))
# Filter contamination rate > 0.8
if (length(which(adipose_Seurat$contamination > 0.8))!=0){
   adipose_Seurat <- subset(adipose_Seurat,contamination <= 0.8)
}
print(paste0("After remove DBL and contamination > 0.8, there are ", nrow(adipose_Seurat@meta.data)," nuclei left"))
individual <- as.data.frame(table(adipose_Seurat@meta.data$donor_id))
write.csv(individual,"individual_distribution_after_doublet_and_contam_removel.csv",row.names=F,quote=F)
# Save the file
write.table(rownames(adipose_Seurat@meta.data),"20.barcode_after_DBL.tsv",row.names = F, quote = F, sep = "\t")
write.table(adipose_Seurat@meta.data,"20.metadata.tsv",row.names = F, quote = F, sep = "\t")

#############################################################################################################
# reprocess the data
options(future.globals.maxSize = 8000 * 1024^2)
adipose_Seurat <- SCTransform(adipose_Seurat,verbose = FALSE,vars.to.regress = c("percent.mt", "percent.ribo","G2M.Score","S.Score"),vst.flavor = "v2")
adipose_Seurat <- Seurat::RunPCA(adipose_Seurat, assay="SCT", verbose = FALSE)
adipose_Seurat <- Seurat::RunUMAP(adipose_Seurat, reduction = "pca",
                       dims = 1:30, verbose = FALSE)
adipose_Seurat <- Seurat::FindNeighbors(adipose_Seurat, reduction = "pca",
                             dims = 1:30, verbose = FALSE)
adipose_Seurat <- Seurat::FindClusters(
   object = adipose_Seurat)
setwd(paste0("./06.seurat_qc/",batch,"/",sample))
saveRDS(adipose_Seurat, "./seurat_rds/Seurat_obj_after_DBL_removal_b4_int_demuxlet.rds")

############################ plot seurat cluster and annotation after QC ################################
setwd("./step04/demuxlet")
# The aim of these two plots is to see whether there is a strong batch effect between individuals before integration
pdf("01.fat_type_distribution_after_remove_DBL.pdf",width=5,height=4)
DimPlot(adipose_Seurat, reduction = "umap",group.by = "Depot",cols = c("VAT" = "#fbb4ae",
                                                                          "SAT" = "#b3cde3"))
dev.off()

pdf("02.sex_type_distribution_after_remove_DBL.pdf",width=5,height=4)
DimPlot(adipose_Seurat, reduction = "umap",group.by = "Gender",cols = c("M" = "#377eb8",
                                                                        "F" = "#e41a1c"))
dev.off()
p2 <- DimPlot(adipose_Seurat, reduction = "umap",group.by = "predicted.celltype.l1",cols = colors)+guides(color=guide_legend(override.aes=list(size=2)))
p3 <- DimPlot(adipose_Seurat, reduction = "umap",group.by = "seurat_clusters",cols = myColors[1:length(unique(adipose_Seurat$seurat_clusters))])+
  guides(color=guide_legend(override.aes=list(size=1)))
pdf("11.annotation_result_after_removing_doublet.pdf",width = 8,height = 10)
p2/p3
dev.off()
# plot individuals
pdf("12.individual_distribution_after_DBL_removel.pdf",width=5,height=4)
DimPlot(adipose_Seurat, reduction = "umap",group.by = "donor_id",cols=individual_col)
dev.off()

# plot individual proportion per celltype

df <- data.frame()
for (i in unique(adipose_Seurat@meta.data$predicted.celltype.l1)) {
  df_tmp <- adipose_Seurat@meta.data[which(adipose_Seurat@meta.data$predicted.celltype.l1==i),c(c("barcode","predicted.celltype.l1","donor_id"))]
  table <- as.data.frame(table(df_tmp$donor_id))
  table$celltype <- i
  table$proportion <- table$Freq/sum(table$Freq)
  df <- rbind(df,table)
}
write.csv(df,"celltype_proportion_and_count_after_DBL_removel.csv",quote=F)
# plot individual proportion per celltype
p <- ggplot(df, aes(x = celltype, y = proportion, fill = Var1)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Proportion of Cell Types by Individual",
       x = "Cell Type",
       y = "Proportion",
       fill = "Individual") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_manual(values=individual_col)
pdf("13.individual_distribution_per_cell_type_after_DBL_removel.pdf",width=5,height=4)
p
dev.off()

df2<-df %>% group_by(celltype) %>% summarise(count=sum(Freq))
df$celltype <- factor(df$celltype,levels=df2$celltype[order(df2$count)])
p <- ggplot(df, aes(x = celltype, y = Freq, fill = Var1)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Count of Cell Types by Individual",
       x = "Cell Type",
       y = "Nuclei count",
       fill = "Individual") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_manual(values=individual_col)
pdf("14.individual_distribution_per_cell_type_by_count_after_DBL_removel.pdf",width=5,height=4)
p
dev.off()

# plot celltype contamination after removing doublets
cont_all <- ggplot(as.data.frame(adipose_Seurat@meta.data),aes(x=contamination))+
        geom_histogram(binwidth = 0.1, fill = "grey", color = "black")+
        ggtitle("Contamination rate of nuclei after doublet removel")+
        theme_classic()
cont_ct <- ggplot(as.data.frame(adipose_Seurat@meta.data),aes(x=contamination))+
        geom_histogram(binwidth = 0.05, fill = "grey", color = "black")+
        facet_wrap(~ predicted.celltype.l1, ncol = 3)+
        ggtitle("Contamination rate of each celltype after doublet removel")+
        theme_classic()
pdf("15.contamination_rate_after_DBL_removal_SCT.pdf",width=5,height=10)
cont_all/cont_ct
dev.off()

# plot celltype proportion after removing doublets
celltype_prop <- as.data.frame(table(adipose_Seurat$predicted.celltype.l1))
celltype_prop$proportion <- celltype_prop$Freq/sum(celltype_prop$Freq)

pdf("16.piechart_of_cellpropotion_after_DBL_remove.pdf",width=5,height=4)
ggplot(celltype_prop, aes(x = "", y = proportion, fill = Var1)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  labs(title = "Cell type proportion") +
  theme_void()+
  scale_fill_manual(values=colors)
dev.off()
individual <- as.data.frame(table(adipose_Seurat@meta.data$donor_id))
#write.csv(individual,"individual_distribution_after_doublet_removel.csv",row.names=F,quote=F)
print("finish plot individual")
print("Doublet part done! Good job😊")


p1 <- DimPlot(adipose_Seurat, reduction = "ref.umap",group.by = "predicted.celltype.l1",cols = myColors[1:length(unique(adipose_Seurat$predicted.celltype.l1))]) + ggtitle("High quality cells")
p2 <- DimPlot(adipose_Seurat, reduction = "umap",group.by = "predicted.celltype.l1",cols = myColors[1:length(unique(adipose_Seurat$predicted.celltype.l1))]) + ggtitle("High quality cells")
pdf("15.annotation_result_after_removing_doublet_and_ambient_SCT.pdf",width = 8,height = 10)
p1/p2
dev.off()
all_ngene <- data.frame(nGene = adipose_Seurat$nFeature_RNA,
                        nUMI = adipose_Seurat$nCount_RNA,
                        pMito = adipose_Seurat$percent.mt)
nGene_tmp <- as.data.frame(adipose_Seurat@meta.data)

nGene_tmp <- nGene_tmp %>%
  filter(!is.infinite(nFeature_RNA), !is.infinite(percent.mt)) %>%
  filter(!is.na(nFeature_RNA), !is.na(percent.mt))
p1 <- ggplot(nGene_tmp, aes(x = nFeature_RNA, y = percent.mt)) +
      geom_point(size = 0.1, color = "grey") +
      stat_density_2d(aes(fill = ..level..), geom = "polygon", color = "white", alpha = 0.5, size = 0.3) +
      scale_fill_viridis_c(option = "plasma", labels = scales::scientific_format()) +
      theme_bw() +
      scale_x_continuous(breaks = seq(0, max(as.numeric(nGene_tmp$nFeature_RNA), na.rm = TRUE), by = 2000), minor_breaks = seq(0, max(as.numeric(nGene_tmp$nFeature_RNA), na.rm = TRUE), by = 500)) +
      scale_y_continuous(breaks = seq(0, max(as.numeric(nGene_tmp$percent.mt), na.rm = TRUE), by = 10)) +
      geom_hline(yintercept = c(0, 5), linetype = "dashed", color = "black") +
      geom_vline(xintercept = 500, linetype = "dashed", color = "red") +
      ggtitle(paste0("n = ", nrow(nGene_tmp)))
p2 <- VlnPlot(adipose_Seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

pdf("17.corrected_distribution_of_NODG_vs_pMito_zoom_in_after_DBL.pdf",width = 8,height = 5)
p1/p2
dev.off()

#################### Step 05: Intergration by each individual ####################
DefaultAssay(adipose_Seurat) <- "RNA"
adipose_Seurat[["RNA"]] <- split(adipose_Seurat[["RNA"]], f = adipose_Seurat$donor_id)
adipose_Seurat <- NormalizeData(adipose_Seurat)
adipose_Seurat <- FindVariableFeatures(adipose_Seurat)
adipose_Seurat <- ScaleData(adipose_Seurat)
adipose_Seurat <- RunPCA(adipose_Seurat)
adipose_Seurat <- FindNeighbors(adipose_Seurat, dims = 1:30, reduction = "pca")
adipose_Seurat<- FindClusters(adipose_Seurat, resolution = 1, cluster.name = "unintegrated_clusters")
adipose_Seurat<- RunUMAP(adipose_Seurat, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
pdf("18.UMAP_before_integration.pdf",width = 8,height=4)
DimPlot(adipose_Seurat, reduction = "umap.unintegrated", group.by = c("predicted.celltype.l1"))
dev.off()
## 03. Perform streamlined (one-line) integrative analysis ###
options(future.globals.maxSize = 3e+09)
adipose_Seurat <- SCTransform(adipose_Seurat,verbose = FALSE,vars.to.regress = c("percent.mt", "percent.ribo","G2M.Score","S.Score"),vst.flavor = "v2")
adipose_Seuratr <- RunPCA(adipose_Seurat, npcs = 30, verbose = F, assay="SCT")
adipose_Seurat <- IntegrateLayers(
  object = adipose_Seurat,
  method = HarmonyIntegration,
  normalization.method = "SCT",
  new.reduction = "harmony",
  verbose = F
)
adipose_Seurat <- FindNeighbors(adipose_Seurat, dims = 1:30, reduction = "harmony")
adipose_Seurat <- FindClusters(adipose_Seurat, resolution = 1,cluster.name = "harmony_clusters")
adipose_Seurat <- RunUMAP(adipose_Seurat, dims = 1:30, reduction = "harmony",reduction.name = "umap.harmony")
adipose_Seurat <- JoinLayers(adipose_Seurat,assay = "RNA")
p1 <- DimPlot(adipose_Seurat, reduction = "ref.umap",group.by = "predicted.celltype.l1",cols = myColors[1:length(unique(adipose_Seurat$predicted.celltype.l1))])
p2 <- DimPlot(adipose_Seurat, reduction = "umap.harmony",group.by = "predicted.celltype.l1",cols = myColors[1:length(unique(adipose_Seurat$predicted.celltype.l1))])
pdf("19.annotation_result_after_removing_doublet.pdf",width = 8,height = 10)
p1/p2
dev.off()

p1 <- DimPlot(
  adipose_Seurat,
  reduction = "umap.harmony",
  group.by = c("predicted.celltype.l1","Gender","donor_id"),
  combine = FALSE, label.size = 2,cols=colors,pt.size = 0.1
)
pdf("21.UMAP_of_harmony_integrative_methods_after_QC_SCT_GT.pdf",width=6,height=4)
p1
dev.off()
print("Integration part done! Good job😊")






           
