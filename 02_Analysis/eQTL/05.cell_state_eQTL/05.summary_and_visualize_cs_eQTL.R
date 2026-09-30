
# Project root is 8 levels up from this script
args0 <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args0[grep("--file=", args0)])
project_root <- normalizePath(file.path(dirname(normalizePath(script_path)),
                                        paste(rep("..", 8), collapse = "/")))

base_data <- file.path(project_root, "data/09.eQTL/03.Saige_QTL/02.output/11.cell_state_eQTL/02.adipogenesis")

sat_pc <- 2
vat_pc <- 3

sat_df <- data.frame()
vat_df <- data.frame()

for (i in 1:22) {
  sat_file <- file.path(base_data, paste0("SAT/output/csQTL/SAT.pPC", sat_pc, ".chr", i, ".diabetic_biased.intersect.tsv"))
  sat_tmp <- read.table(sat_file, sep = "\t", header = TRUE)
  sat_df <- rbind(sat_df, sat_tmp)

  vat_file <- file.path(base_data, paste0("VAT/output/csQTL/VAT.pPC", vat_pc, ".chr", i, ".diabetic_biased.intersect.tsv"))
  vat_tmp <- read.table(vat_file, sep = "\t", header = TRUE)
  vat_df <- rbind(vat_df, vat_tmp)
}

sat_df$FDR <- p.adjust(sat_df$pvalue,method = "BH")
vat_df$FDR <- p.adjust(vat_df$pvalue,method = "BH")

sat_sig <- sat_df[which(sat_df$FDR <= 0.05),]
sat_sig$eQTL <- paste0(sat_sig$gene_id,"-",sat_sig$snp_lead)
vat_sig <- vat_df[which(vat_df$FDR <= 0.05),]
vat_sig$eQTL <- paste0(vat_sig$gene_id,vat_sig$snp_lead)

merge <- rbind(sat_sig,vat_sig)
write.table(merge, file.path(base_data, "SAT_pc2_VAT_pc3.sig_cseQTL.ALL.FDR0.05.tsv"), sep = "\t", quote= F, row.names = F)

shared <- merge(sat_sig, vat_sig, by = c("gene_id", "snp_id"))
write.table(shared, file.path(base_data, "SAT_pc2_VAT_pc3.sig_cseQTL.shared.FDR0.05.tsv"), sep = "\t", quote= F, row.names = F)





