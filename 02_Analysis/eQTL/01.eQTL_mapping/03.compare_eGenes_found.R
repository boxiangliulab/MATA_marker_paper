# Project root is 8 levels up from this script
args0 <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args0[grep("--file=", args0)])
project_root <- normalizePath(file.path(dirname(normalizePath(script_path)),
                                        paste(rep("..", 8), collapse = "/")))

filelist <- list.files("/home/users/nus/e1101943/scratch/project/00.adipose/data/09.eQTL/01.SAIGE_QTL/03.formal_calling/SAT/")
dat <- data.frame()
for (i in filelist){
        setwd(file.path(project_root, "data/09.eQTL/03.Saige_QTL/02.output/03.formal_calling/ancestry_specific_version/SAT/01.egene_list"))
        folder_path = file.path(project_root, "data/09.eQTL/03.Saige_QTL/02.output/03.formal_calling/ancestry_specific_version/SAT", i)
        matched_files <- list.files(path = folder_path,
                           pattern = "genePval",
                           full.names = TRUE)
        if (length(matched_files) != 0){
        df <- data.frame()
        for (j in matched_files){
                df_tmp <- read.table(j,header = T)
                df_tmp$depot <- "SAT"
                df_tmp$celltype <- i
                df <- rbind(df,df_tmp)
                }
        df$FDR <- p.adjust(df$ACAT_p, method = "BH")
        write.table(df,paste0("SAT.",i,".gene_all_list.tsv"),sep = "\t",quote=F,row.names = F)
        write.table(df[which(df$FDR <= 0.05),],paste0("SAT.",i,".egene_list.tsv"),sep = "\t",quote=F,row.names = F)
        ngene = length(which(df$FDR <= 0.05))
        print(paste0(i," found ",ngene, " egenes"))
        dat_tmp <- data.frame(celltype = i,n_eGene = ngene)
        dat <- rbind(dat,dat_tmp)
        print(paste0(i," done!"))
        }
}
