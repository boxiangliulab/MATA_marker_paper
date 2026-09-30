### SATA  Saige eQTL analysis:
# version 1: 26 December 2024
# version 2: 1st April 2025. Modify: Input from raw count data from scanpy, create obj from raw count
### Function: Generate phenotype for saigeQTL
### Author: Yihan Tong
library(Seurat)
library(tidyverse)
library(stringr)
library(cowplot)
library(patchwork)
library(fastDummies)
library(sctransform)
library(getopt)
library(data.table)

command=matrix(c("depot","d",1,"character",
                 "celltype","c",1,"character",
                 "help","h",0,"logical"),byrow=T,ncol=4)
args=getopt(command)
if( !is.null(args$help) || is.null(args$celltype) || is.null(args$depot))
{ cat(paste(getopt(command,usage=T),"\n"))
  q() }

celltype = args$celltype
depot = args$depot

# Project root is 8 levels up from this script
args0 <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args0[grep("--file=", args0)])
project_root <- normalizePath(file.path(dirname(normalizePath(script_path)),
                                        paste(rep("..", 8), collapse = "/")))

sat_path = read.table("/home/users/nus/e1101943/scratch/project/00.adipose/data/07.full_scale_snRNA/08.integration/all/scvi/03.count_matrix/SAT/01.raw_count/00.SAT_celltype_count.path")
vat_path = read.table("/home/users/nus/e1101943/scratch/project/00.adipose/data/07.full_scale_snRNA/08.integration/all/scvi/03.count_matrix/VAT/01.raw_count/00.VAT_celltype_count.path")
metadata = read.table(file.path(project_root, "data/08.integration_snRNA_object/01.freeze_version/metadata_after_libraryQC_after_ctQC_freeze.tsv"), sep = "\t", header = T)
metadata <- metadata[-which(metadata$Race == "Caucasian"),]

## 2025.4.14 version: Filter out relateness indivdivual (Filter the less nuclei count individual in the relateness pair)
relateness_pair <- read.csv("/home/project/11003054/e1353486/project/adipose/fullgenotype/new_version/08.relatedness/close_relative.csv",header = T)
relate_to_remove <- relateness_pair %>%
  #SAT
  mutate(
    Depot = "SAT",
    donor_id = ifelse(IID1_SAT_cell_number < IID2_SAT_cell_number, IID1, IID2)
  ) %>%
  select(Depot, donor_id) %>%
  #VAT
  bind_rows(
    relateness_pair %>%
      mutate(
        Depot = "VAT",
        donor_id = ifelse(IID1_VAT_cell_number < IID2_VAT_cell_number, IID1, IID2)
      ) %>%
      select(Depot, donor_id)
  )

# Filtered relateness
metadata_sat <- metadata[which(metadata$Depot == "SAT"),]
relate_to_remove_sat <- relate_to_remove[which(relate_to_remove$Depot == "SAT"),]
metadata_sat <- metadata_sat[-which(metadata_sat$donor_id %in% relate_to_remove_sat$donor_id),]

metadata_vat <- metadata[which(metadata$Depot == "VAT"),]
relate_to_remove_vat <- relate_to_remove[which(relate_to_remove$Depot == "VAT"),]
metadata_vat <- metadata_vat[-which(metadata_vat$donor_id %in% relate_to_remove_vat$donor_id),]
metadata <- rbind(metadata_sat,metadata_vat)
print(paste0("There are ",nrow(metadata),"in metadata", "unique line number is:", length(unique(metadata$barcode))))
## subset donors with => 10 nuclei in this cell type
df3 <- metadata %>%
  group_by(Depot, donor_id,Annotation_level2) %>%
  summarise(count = n(), .groups = "drop")
df3.filter <- df3[which(df3$count >=10),]
df3.sat <- df3[which(df3$Depot == "SAT"),]
df3.vat <- df3[which(df3$Depot == "VAT"),]

df3.sat <- df3.sat %>%
  filter(count >= 10) %>%
  group_by(Annotation_level2) %>%
  summarise(n_donors = n_distinct(donor_id), .groups = "drop")

df3.vat <- df3.vat %>%
  filter(count >= 10) %>%
  group_by(Annotation_level2) %>%
  summarise(n_donors = n_distinct(donor_id), .groups = "drop")

## Select celltypes with >= 10 donors
sat_celltype <- df3.sat$Annotation_level2[which(df3.sat$n_donors >= 10)]
print(paste0("There are ", length(sat_celltype)," cell types tested in SAT"))

## Check if all cell types are in the metadata
for (i in sat_celltype) {
        if (i %in% metadata$Annotation_level2){
                next
                } else {
                        message(paste("SAT: Cell type", i, "not found in metadata$Annotation_level2"))
                        }
        }

vat_celltype <- df3.vat$Annotation_level2[which(df3.vat$n_donors >= 10)]
print(paste0("There are ", length(vat_celltype)," cell types tested in VAT"))

for (i in vat_celltype) {
        if (i %in% metadata$Annotation_level2){
                next
                } else {
                        message(paste("VAT: Cell type", i, "not found in metadata$Annotation_level2"))
                        }
        }
get_saige_phenotype <- function(depot,celltype,n_gpc,n_ppc) {
    print(paste0("[] Start processing ", celltype))
    options(future.globals.maxSize = 1000000 * 1024^2)
    #count <- read.table(paste0("/home/users/nus/e1101943/scratch/project/00.adipose/data/07.full_scale_snRNA/08.integration/all/scvi/03.count_matrix/",depot,"/01.raw_count/",celltype,".raw_count.tsv"),sep = "\t",header = T,row.names= 1,check.names=F)
    file_path <- paste0("/home/users/nus/e1101943/scratch/project/00.adipose/data/07.full_scale_snRNA/08.integration/all/scvi/03.count_matrix/",
                     depot, "/01.raw_count/", celltype, ".raw_count.tsv")
    count <- fread(file_path, sep = "\t") |>
        as.data.frame()
    rownames(count) <- count[,1]
    count <- count[,-1]
    #count <- as.data.frame(count)
    count <- t(count) %>% as.matrix()
    count <- count[,which(colnames(count) %in% metadata$barcode)]
    obj.filter <- CreateSeuratObject(count,project = celltype)
    print(paste0("There are ", nrow(obj.filter@meta.data), " nuclei in ", celltype))
    print(obj.filter)
    meta_sub <- metadata[which(metadata$barcode %in% rownames(obj.filter@meta.data)),]

    rownames(meta_sub) <- meta_sub$barcode
    obj.filter <- AddMetaData(obj.filter,meta_sub)
    print(head(obj.filter))
    setwd(file.path(project_root, "data/09.eQTL/03.Saige_QTL/01.input/01.phenotype", depot))
    #setwd(paste0("/home/users/nus/e1101943/scratch/project/00.adipose/data/09.eQTL/01.SAIGE_QTL/01.input/01.phenotype/",depot))
    if(!dir.exists(celltype)){dir.create(celltype)}
    setwd(celltype)
    #Normalize the data
    obj.filter <- NormalizeData(object = obj.filter, normalization.method = "LogNormalize", scale.factor = 10000)
    #Find HVGs
    obj.filter <- FindVariableFeatures(obj.filter, selection.method = "vst", nfeatures = 2000)#Find Clusters
    #Center and scale the data
    all.genes <- rownames(obj.filter)
    obj.filter <- ScaleData(object = obj.filter,features = all.genes )
    #Run PCA on the data
    obj.filter <- RunPCA(object = obj.filter, features = VariableFeatures(object = obj.filter), verbose = FALSE)
    #Find Neighbors
    obj.filter <- FindNeighbors(obj.filter, reduction = "pca", dims = 1:20)
    #Find obj.filter
    obj.filter <- FindClusters(obj.filter, resolution = 1)
    # Filter donors has >= 10 nuclei in this cell type
    df_donors_of_interest <- as.data.frame(table(obj.filter[["donor_id"]]))
    vec_donors_of_interest <- as.character(df_donors_of_interest$donor_id[which(df_donors_of_interest$Freq >= 10)])
    # keep genes expressed more than 2% of nuclei
    matrix_SATA_Asian <- GetAssayData(object = obj.filter,assay = "RNA", layer = "counts")
    vec_genes_numberofcells <- apply(matrix_SATA_Asian, MARGIN = 1, function (x) length(which(x > 0)))
    vec_genes_to_keep <- vec_genes_numberofcells[vec_genes_numberofcells > dim(matrix_SATA_Asian)[2]*0.02]
    vec_genes_to_keep <- names(vec_genes_to_keep)

    eQTL_subset <- subset(obj.filter,donor_id %in% vec_donors_of_interest)
    eQTL_matrix <- GetAssayData(object = eQTL_subset,layer= "counts",assay="RNA") %>% as.data.frame() %>% t()
    vec_genes_to_keep <- vec_genes_to_keep[which(vec_genes_to_keep %in% colnames(eQTL_matrix))]
    eQTL_matrix <- eQTL_matrix[,vec_genes_to_keep] %>% as.data.frame()
    eQTL_matrix$barcode <- rownames(eQTL_matrix)

     # Wrangle the covariants: sex, age, BMI, log10(total UMI+1), genotype PC (PC), diabetic
    eQTL_subset$log1p_UMI <- log10(eQTL_subset$nCount_RNA + 1)
    eQTL_subset@meta.data$sex <- 0
    eQTL_subset@meta.data$sex <- ifelse(eQTL_subset@meta.data$Gender=="M",0,1)
    eQTL_subset@meta.data$diabetic <- 0
    eQTL_subset@meta.data$diabetic <- ifelse(eQTL_subset@meta.data$Self_defined_DM=="T2DM",1,0)

    # read genotype pc in
    # genotype_pc <- read.csv("/home/project/11003054/e1353486/project/adipose/fullgenotype/PCA/300_sata_pca_ethnicity_result.csv") v1
    genotype_pc <-  read.csv("/home/project/11003054/e1353486/project/adipose/fullgenotype/Ancestry_specific/01.PCA/296_sata_pca_ethnicity_result.csv") # ancestry specific version
    rownames(genotype_pc) <- genotype_pc$FID

    if (length(which(genotype_pc$FID %in% vec_donors_of_interest)) == length(vec_donors_of_interest)){
            gPC_colum <- c("FID",paste0("PC",1:n_gpc))
            genotype_pc.filter <- genotype_pc[which(genotype_pc$FID %in% vec_donors_of_interest),gPC_colum]
            } else {
                    message("Oops, there is something wrong with genotype name")
                    }
     gPC_colum <- paste0("PC",1:n_gpc)
    column_filter <- c("barcode","donor_id","sex","diabetic","Age","BMI","log1p_UMI","percent.mt",gPC_colum)
    # get matadata
    meta.filter <- as.data.frame(eQTL_subset@meta.data)
    meta.filter <- merge(meta.filter,genotype_pc.filter,by.x = "donor_id",by.y = "FID")
    meta.filter <- meta.filter[,column_filter]

    # get gene expression PC for each nuclei
    ppc_column <- c(paste0("PC_",1:n_ppc))
    pPC <- FetchData(obj.filter,ppc_column)
    pPC$cell <- row.names(pPC)
    row.names(pPC) <- NULL
    meta.filter <- merge(meta.filter,pPC,by.x = "barcode",by.y = "cell")
    colnames(meta.filter) <- gsub("^PC_(\\d+)", "pPC\\1",
                              gsub("^PC(\\d+)", "gPC\\1", colnames(meta.filter)))# change PC1, PC2, ... to gPC1, gPC2, ...,# change PC_1, PC_2, ... to pPC1, pPC2, ...
    # merge
    merge <- merge(meta.filter, eQTL_matrix, by = "barcode")
    merge$barcode <- NULL
    #write.table(merge,paste0(celltype,"_saige_phenotype_sct_raw_gene_name.tsv"),sep = "\t",quote = F, row.names=F )
    names(merge) <- gsub("-","_",names(merge))
    write.table(merge,paste0(celltype,"_saige_phenotype.RNA.safe_gene_name.","gPC",n_gpc,"_pPC",n_ppc,".tsv"),sep = "\t",quote = F, row.names=F )
    #write.table(merge,paste0(celltype,"_saige_phenotype.raw.safe_gene_name.tsv"),sep = "\t",quote = F, row.names=F )
    # in case saige cannot recognize
    #sed '1s/\./\\\\./g'  Adipocyte_saige_phenotype_raw_gene_name.tsv > Adipocyte_saige_phenotype_safe_gene_name.tsv
    print(paste0("😊 Finish ", celltype))
}

### SAT
library(purrr)
library(furrr)
library(future)
plan(multisession)

# SAT
#for (celltype in sat_celltype){
#                get_saige_phenotype(depot = "SAT",celltype = celltype,n_gpc=8,n_ppc=10)
#        }

## VAT
#for (celltype in vat_celltype){
#                get_saige_phenotype(depot = "VAT",celltype = celltype,n_ppc=10,n_gpc=6)
#}

if (depot=="SAT"){
        get_saige_phenotype(depot = "SAT",celltype = celltype,n_gpc=8,n_ppc=10)
        } else if (depot=="VAT") {
                get_saige_phenotype(depot = "VAT",celltype = celltype,n_ppc=10,n_gpc=6)
                }






