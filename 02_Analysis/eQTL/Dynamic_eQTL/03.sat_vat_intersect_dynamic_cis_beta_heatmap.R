library(data.table)
library(dplyr)
library(tidyr)
library(pheatmap)


# 1. Read dynamic eQTL results and calculate FDR

sig_dynamic_eqtl <- function(depot){
  input <- paste0(
    "/project/adipose/dynamic_qtl/02.output/",
    depot,
    "/dynamic_eQTL/"
  )
  
  filelist <- list.files(input, pattern = "\\.tsv$", full.names = TRUE)
  
  df <- rbindlist(lapply(filelist, fread), fill = TRUE)
  df[, p_dynamic_fdr := p.adjust(p_dynamic, method = "BH")]
  df[, depot := depot]
  
  return(df)
}

sat <- sig_dynamic_eqtl("SAT")
vat <- sig_dynamic_eqtl("VAT")

merge <- rbindlist(list(sat, vat), fill = TRUE)
merge[, lead_eqtl := paste0(gene_id, "-", snp_lead)]
merge_sig <- merge[p_dynamic_fdr <= 0.05]

# 只保留 lead SNP 行
merge_sig <- merge_sig %>%
  filter(snp_id == snp_lead)
top_df <- merge_sig %>%
  as.data.frame() %>%
  distinct(depot, gene_id, snp_lead, .keep_all = TRUE) %>%
  group_by(gene_id, snp_lead) %>%
  filter(n_distinct(depot) == 2) %>%
  ungroup()
top_df <- top_df %>%
  distinct(gene_id, snp_lead, .keep_all = TRUE)

# =========================================================
# 2. Take top100 unique gene-snp signals
# =========================================================
top100 <- top_df
#<- top_df %>%
#  as.data.frame() %>%
#  arrange(p_dynamic_fdr, p_dynamic) %>%
#  distinct(gene_id, snp_lead, .keep_all = TRUE) %>%
#  slice_head(n = 100)

top100$row_id <- paste0(top100$gene_id, "_", top100$snp_lead)


# =========================================================
# 2. column order
# =========================================================
celltype_order <- c(
  "Early_preadipocyte",
  "Committed_preadipocyte",
  "Adipocyte"
)

col_order <- c(
  "SAT_Early_preadipocyte",
  "SAT_Committed_preadipocyte",
  "SAT_Adipocyte",
  "VAT_Early_preadipocyte",
  "VAT_Committed_preadipocyte",
  "VAT_Adipocyte"
)

# =========================================================
# 3. function: get raw BETA from cis file
# =========================================================
get_beta <- function(gene, snp, depot, celltype){
  f <- paste0(
    "/project/01.adipose_sqtl/data/09.eQTL/03.Saige_QTL/02.output/03.formal_calling/ancestry_specific_version/",
    depot, "/", celltype, "/", gene, "_cis"
  )
  
  x <- fread(f)
  hit <- x[MarkerID == snp]
  hit$BETA[1]
}

# =========================================================
# 4. build beta matrix
# =========================================================
beta_mat <- matrix(
  NA_real_,
  nrow = nrow(top100),
  ncol = length(col_order)
)

rownames(beta_mat) <- top100$row_id
colnames(beta_mat) <- col_order

get_beta <- function(gene, snp, depot, celltype){
  f <- paste0(
    "/project/01.adipose_sqtl/data/09.eQTL/03.Saige_QTL/02.output/03.formal_calling/ancestry_specific_version/",
    depot, "/", celltype, "/", gene, "_cis"
  )
  
  if (!file.exists(f)) return(NA_real_)
  
  x <- fread(f)
  hit <- x[MarkerID == snp]
  
  if (nrow(hit) == 0) return(NA_real_)
  
  return(as.numeric(hit$BETA[1]))
}
for (i in seq_len(nrow(top100))) {
  
  g <- top100$gene_id[i]
  s <- top100$snp_lead[i]
  
  beta_mat[i, "SAT_Early_preadipocyte"] <- get_beta(g, s, "SAT", "Early_preadipocyte")
  beta_mat[i, "SAT_Committed_preadipocyte"] <- get_beta(g, s, "SAT", "Committed_preadipocyte")
  beta_mat[i, "SAT_Adipocyte"] <- get_beta(g, s, "SAT", "Adipocyte")
  
  beta_mat[i, "VAT_Early_preadipocyte"] <- get_beta(g, s, "VAT", "Early_preadipocyte")
  beta_mat[i, "VAT_Committed_preadipocyte"] <- get_beta(g, s, "VAT", "Committed_preadipocyte")
  beta_mat[i, "VAT_Adipocyte"] <- get_beta(g, s, "VAT", "Adipocyte")
}
# =========================================================
# 5. save matrix
# =========================================================
write.table(
  cbind(row_id = rownames(beta_mat), as.data.frame(beta_mat)),
  "/project/adipose/dynamic_qtl/02.output/dynamic_cis_beta_matrix.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# =========================================================
# 6. plot heatmap
# =========================================================
max_abs <- max(abs(beta_mat), na.rm = TRUE)
breaks <- seq(-max_abs, max_abs, length.out = 101)
cols <- colorRampPalette(c("#3B4CC0", "white", "#E31A1C"))(100)

beta_mat_plot <- beta_mat[, col_order, drop = FALSE]
beta_mat_clust <- beta_mat_plot
beta_mat_clust[is.na(beta_mat_clust)] <- rowMeans(beta_mat_clust, na.rm = TRUE)[row(beta_mat_clust)][is.na(beta_mat_clust)]
row_hc <- hclust(dist(beta_mat_clust))

pdf("/project/adipose/dynamic_qtl/02.output/sat_vat_intersect_dynamic_cis_beta_heatmap.pdf", width = 7, height = 30)
pheatmap(
  beta_mat_plot,
  cluster_rows = row_hc,
  cluster_cols = FALSE,
  color = cols,
  breaks = breaks,
  gaps_col = 3,
  na_col = "grey90",
  border_color = NA,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 5,
  fontsize_col = 10
)
dev.off()

