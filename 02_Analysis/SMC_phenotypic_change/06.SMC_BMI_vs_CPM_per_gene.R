suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(Seurat)
  library(edgeR)
  library(ggplot2)
  library(ggrepel)
  library(Matrix)
  library(getopt)
})

# =======================
# Args
# =======================
command <- matrix(c(
  "celltype", "c", 1, "character",
  "help",     "h", 0, "logical"
), byrow = TRUE, ncol = 4)

args <- getopt(command)

if (!is.null(args$help) || is.null(args$celltype)) {
  cat(paste(getopt(command, usage = TRUE), "\n"))
  q(status = 1)
}

celltype <- args$celltype

# Project root is 7 levels up from this script
args0 <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args0[grep("--file=", args0)])
project_root <- normalizePath(file.path(dirname(normalizePath(script_path)),
                                        paste(rep("..", 7), collapse = "/")))

# =======================
# Config
# =======================
min_cells_per_donor <- 10

cols <- c("SAT" = "#BBCFC3", "VAT" = "#FBCE9E")

# =======================
# Paths
# =======================
metadata_path <- file.path(project_root, "data/08.integration_snRNA_object/01.freeze_version/metadata_after_libraryQC_after_ctQC_freeze.tsv")

base_count_dir <- "/home/users/nus/e1101943/scratch/project/00.adipose/data/07.full_scale_snRNA/08.integration/all/scvi/03.count_matrix"
count_path_sat <- file.path(base_count_dir, "SAT", "01.raw_count", paste0(celltype, ".raw_count.tsv"))
count_path_vat <- file.path(base_count_dir, "VAT", "01.raw_count", paste0(celltype, ".raw_count.tsv"))

sig_table_dir <- file.path(project_root, "data/21.differential_expressed_genes/01.pseudobulk/04.only_VAT_vs_SAT_pseudobulk/linear_model/interation_term/tables")
sig_file <- file.path(sig_table_dir, paste0(celltype, ".DepotBMI.sig.limma.tsv"))

plot_dir <- file.path(project_root,
  "data/21.differential_expressed_genes/01.pseudobulk/04.only_VAT_vs_SAT_pseudobulk/linear_model/interation_term/plots",
  celltype
)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# =======================
# Sanity check
# =======================
message("[Args] celltype=", celltype)
if (!file.exists(count_path_sat)) stop("SAT count not found: ", count_path_sat)
if (!file.exists(count_path_vat)) stop("VAT count not found: ", count_path_vat)
if (!file.exists(sig_file)) stop("Sig table not found: ", sig_file)

# =======================
# Read sig genes
# =======================
sig_df <- fread(sig_file, sep = "\t", header = TRUE, data.table = FALSE)
sig_genes <- unique(sig_df$gene)
message("[SigGenes] ", length(sig_genes), " genes to plot for ", celltype)

if (length(sig_genes) == 0) {
  message("[!] No sig genes. Exiting.")
  q(status = 0)
}

# =======================
# Load metadata
# =======================
metadata <- fread(metadata_path, sep = "\t", header = TRUE, data.table = FALSE)

# =======================
# Helper: read raw_count.tsv -> Seurat
# =======================
read_count_as_seurat <- function(fp) {
  m <- fread(fp, sep = "\t", header = TRUE, check.names = FALSE, data.table = FALSE)
  first_col <- colnames(m)[1]
  rownames(m) <- m[[first_col]]
  m[[first_col]] <- NULL
  m <- t(as.matrix(m))
  CreateSeuratObject(counts = m)
}

# =======================
# Read and merge
# =======================
message("[Read] SAT count: ", count_path_sat)
seu_sat <- read_count_as_seurat(count_path_sat)
message("[Read] VAT count: ", count_path_vat)
seu_vat <- read_count_as_seurat(count_path_vat)
seu <- merge(seu_sat, y = seu_vat)

# Add metadata
meta_sub <- metadata[metadata$barcode %in% rownames(seu@meta.data), ]
rownames(meta_sub) <- meta_sub$barcode
seu <- AddMetaData(seu, meta_sub)
seu <- JoinLayers(seu)

# Filter to Asian donors with both depots
seu <- subset(
  seu,
  Race %in% c("Malay", "Chinese", "Indian") &
    Depot %in% c("SAT", "VAT")
)

if (ncol(seu) == 0) {
  message("[!] No cells after filtering. Exiting.")
  q(status = 0)
}

seu$Depot <- factor(seu$Depot, levels = c("SAT", "VAT"))

# =======================
# Donor filtering: min cells in each depot, paired donors
# =======================
df_cnt <- seu@meta.data %>%
  count(donor_id, Depot, name = "n")

donors_sat_ok <- df_cnt %>% filter(Depot == "SAT", n >= min_cells_per_donor) %>% pull(donor_id)
donors_vat_ok <- df_cnt %>% filter(Depot == "VAT", n >= min_cells_per_donor) %>% pull(donor_id)
donors_keep <- intersect(donors_sat_ok, donors_vat_ok)

message("[DonorFilter] paired donors kept: ", length(donors_keep))

if (length(donors_keep) < 3) {
  message("[!] Too few donors. Exiting.")
  q(status = 0)
}

seu <- subset(seu, donor_id %in% donors_keep)

# =======================
# Pseudobulk per donor x Depot
# =======================
pb_counts <- AggregateExpression(
  seu,
  assays = "RNA",
  group.by = c("donor_id", "Depot"),
  slot = "counts"
)$RNA
df_counts <- as.data.frame(pb_counts)

# =======================
# Pseudobulk metadata
# =======================
meta_pb <- seu[[]] %>%
  dplyr::select(donor_id, Depot, BMI) %>%
  dplyr::distinct(donor_id, Depot, .keep_all = TRUE)

meta_pb$pb_id <- paste(meta_pb$donor_id, meta_pb$Depot, sep = "_")
pb_ids <- colnames(df_counts)
meta_pb <- meta_pb[match(pb_ids, meta_pb$pb_id), ]
stopifnot(all(meta_pb$pb_id == pb_ids))

# =======================
# CPM and log2(CPM+1)
# =======================
dge <- DGEList(counts = df_counts)
dge <- calcNormFactors(dge)
cpm_mat <- cpm(dge, normalized.lib.sizes = TRUE, log = FALSE)
log2cpm <- log2(cpm_mat + 1)

# =======================
# Plot each sig gene: BMI vs log2(CPM+1)
# =======================
sig_genes_avail <- intersect(sig_genes, rownames(log2cpm))
message("[Plot] ", length(sig_genes_avail), " sig genes found in expression matrix")

if (length(sig_genes_avail) == 0) {
  message("[!] No sig genes in expression matrix. Exiting.")
  q(status = 0)
}

for (g in sig_genes_avail) {
  expr_vec <- log2cpm[g, ]

  plot_df <- data.frame(
    pb_id = names(expr_vec),
    Expression = as.numeric(expr_vec),
    BMI = meta_pb$BMI,
    Depot = as.character(meta_pb$Depot),
    donor_id = meta_pb$donor_id,
    stringsAsFactors = FALSE
  )
  plot_df$Depot <- factor(plot_df$Depot, levels = c("SAT", "VAT"))

  # Get sig stats for subtitle
  gene_row <- sig_df[sig_df$gene == g, , drop = FALSE]
  if (nrow(gene_row) > 0) {
    beta_txt <- sprintf("%.4f", gene_row$beta[1])
    fdr_txt  <- formatC(gene_row$fdr[1], format = "e", digits = 2)
    subtitle_txt <- paste0("DepotVAT:BMI beta=", beta_txt, "  FDR=", fdr_txt)
  } else {
    subtitle_txt <- ""
  }

  p <- ggplot(plot_df, aes(x = BMI, y = Expression, color = Depot)) +
    geom_point(size = 2, alpha = 0.8) +
    geom_smooth(method = "lm", se = TRUE, alpha = 0.2, linewidth = 0.8) +
    scale_color_manual(values = cols) +
    labs(
      title = paste0(celltype, ": ", g),
      subtitle = subtitle_txt,
      x = "BMI",
      y = "log2(CPM+1)"
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.subtitle = element_text(size = 8),
      legend.position = "top"
    )

  pdf(file.path(plot_dir, paste0(g, ".BMI_vs_log2CPM.pdf")), width = 5, height = 4)
  print(p)
  dev.off()
}

message("[Done] ", celltype, " | ", length(sig_genes_avail), " plots saved to ", plot_dir)
