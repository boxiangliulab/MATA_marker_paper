library(data.table)
library(tidyverse)

# Project root is 8 levels up from this script
args0        <- commandArgs(trailingOnly = FALSE)
script_path  <- sub("--file=", "", args0[grep("--file=", args0)])
project_root <- normalizePath(file.path(dirname(normalizePath(script_path)),
                                        paste(rep("..", 8), collapse = "/")))

bed <- read.table(file.path(project_root, "data/09.eQTL/02.saige_vs_QTLtools/gene_positions_safe_name.bed"), sep = "\t", header = T)
#lead_eqtl <- read.table("/home/users/nus/e1101943/scratch/project/00.adipose/data/09.eQTL/01.SAIGE_QTL/05.snpEff/SAT/02.egene_per_celltype.tsv",sep = "\t",header = T)
#position <- read.table("/home/users/nus/e1101943/scratch/project/00.adipose/data/09.eQTL/01.SAIGE_QTL/05.snpEff/SAT/lead_eSNP_position.tsv",sep = " ",header = T)
data <- read.table(file.path(project_root, "data/09.eQTL/03.Saige_QTL/02.output/03.formal_calling/ancestry_specific_version/All_eQTL_SAT_and_VAT_with_beta.tsv"), sep = "\t", header=T)
# Convert 0-based bed into 1-based to match with the vcf coordinate
bed$start <-  bed$start + 1
bed$gene_name <- gsub("_","-", bed$gene_name)
# Merge data frames to match the coordiante of rsID
#data <- merge(lead_eqtl,position,by.x = "top_MarkerID",by.y = "ID")
data$CHROM <- str_split_fixed(data$MarkerID,":",4)[,1]
data$POS <- str_split_fixed(data$MarkerID,":",4)[,2]
data$POS <- as.numeric(data$POS)
data$position <- paste0(data$CHROM,":",data$POS)
#data$eqtl <- paste0(data$gene,"-",data$position)
#data$position <- sub("^chr([0-9XYM]+):([0-9]+):.*$", "chr\\1:\\2", data$MarkerID)
data <- data[,c("gene","position","eQTL2","CHROM","POS")] %>% unique()
names(data) <- c("gene","position","eqtl","CHROM","POS")
df <- data.frame()
for (i in 1:nrow(data)){
        gene = data$gene[i]
        var_pos = data$POS[i]
        df_tmp <- bed[which(bed$gene_name == gene),]
        if (nrow(df_tmp)==1){
                df_tmp = df_tmp
                } else {
                        df_tmp = df_tmp[1,]
                        }
        if (df_tmp$strand == "+") {
                distance_tss = var_pos - df_tmp$start # distance to tss
                distance_tes = var_pos - df_tmp$end # distance to tes
                eqtl = data$eqtl[i]
                chr = data$CHROM[i]
                pos = data$POS[i]
                tmp <- data.frame(distance_tss = distance_tss, distance_tes = distance_tes,eqtl = eqtl,gene = gene, chr = chr, pos = pos)
                df <- rbind(df,tmp)
                } else if (df_tmp$strand == "-") {
                        distance_tss = df_tmp$end - var_pos # distance to tss
                        distance_tes = df_tmp$start - var_pos # distance to tes
                        eqtl = data$eqtl[i]
                        chr = data$CHROM[i]
                        pos = data$POS[i]
                        tmp <- data.frame(distance_tss = distance_tss, distance_tes = distance_tes, eqtl = eqtl,gene = gene, chr = chr, pos = pos)
                        df <- rbind(df,tmp)
                        }
       
                
                }
#write.table(df,"/home/users/nus/e1101943/scratch/project/00.adipose/data/09.eQTL/01.SAIGE_QTL/06.tss_region/leadSNP_distance_to_TSS_and_TES_all_celltype.tsv",sep = "\t",quote = F, row.names = F)
write.table(df, file.path(project_root, "data/09.eQTL/03.Saige_QTL/02.output/05.tss_region/ALL_leadSNP_distance_to_TSS_and_TES_all_celltype.tsv"), sep = "\t", quote = F, row.names = F)








