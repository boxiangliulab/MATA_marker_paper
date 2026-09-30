# =============================================================
# Intersection pathways (stress/immune) + UMAP expression
# SAT: PC2, VAT: PC3
# =============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(Matrix)
  library(edgeR)
})

# Project root is 8 levels up from this script
args0 <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args0[grep("--file=", args0)])
project_root <- normalizePath(file.path(dirname(normalizePath(script_path)),
                                        paste(rep("..", 8), collapse = "/")))

# ---------- paths ----------
base_data <- file.path(project_root, "data/09.eQTL/03.Saige_QTL/02.output/11.cell_state_eQTL/02.adipogenesis")
out_dir   <- file.path(base_data, "UMAP_stress_immune_plots")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

sat_pc <- 2
vat_pc <- 3

expr_pal <- c("#053061","#2166ac","#4393c3","#92c5de","#d1e5f0",
              "#f7f7f7","#fddbc7","#f4a582","#d6604d","#b2182b","#67001f")

# ---------- 1. load GSEA ----------
message("Loading GSEA files...")

load_gsea <- function(depot, pc) {
  f <- file.path(base_data, depot, "output", "GSEA",
                 sprintf("adipogenesis.%s.PC%d.gseGO_BP.csv", depot, pc))
  read_csv(f, show_col_types = FALSE) %>% mutate(depot = depot)
}

sat_df <- load_gsea("SAT", sat_pc)
vat_df <- load_gsea("VAT", vat_pc)

# ---------- 2. intersection pathways (same logic as ridge plot script) ----------
sig_cut <- 0.05

gsea <- bind_rows(sat_df, vat_df) %>%
  transmute(
    depot,
    Description     = as.character(Description),
    NES             = as.numeric(NES),
    padj            = as.numeric(p.adjust),
    core_enrichment
  ) %>%
  filter(is.finite(NES), is.finite(padj))

shared_sig_terms <- gsea %>%
  filter(padj <= sig_cut) %>%
  distinct(Description, depot) %>%
  count(Description) %>%
  filter(n == 2) %>%
  pull(Description)

message(sprintf("Shared significant terms: %d", length(shared_sig_terms)))

gsea_shared <- gsea %>% filter(Description %in% shared_sig_terms)

# ---------- 3. extract stress & immune genes from intersection ----------
message("Extracting stress and immune pathway genes...")

target_keywords <- c("stress", "immune")

# genes from intersection pathways matching the keyword
target_genes_combined <- map(target_keywords, function(kw) {
  gsea_shared %>%
    filter(str_detect(tolower(Description), kw)) %>%
    pull(core_enrichment) %>%
    str_split("/") %>%
    unlist() %>%
    trimws() %>%
    unique()
})
names(target_genes_combined) <- target_keywords

message(sprintf("Stress genes: %d", length(target_genes_combined$stress)))
message(sprintf("Immune genes: %d", length(target_genes_combined$immune)))

# ---------- 4. load count matrix + log(CPM+1) ----------
message("Loading count matrices...")

load_counts_logcpm <- function(depot) {
  base <- file.path(base_data, depot)
  mat  <- readMM(file.path(base, "counts.mtx"))      # genes x cells
  barcodes <- read.table(file.path(base, "barcodes.tsv"),
                         header = FALSE, stringsAsFactors = FALSE)$V1
  genes    <- read.table(file.path(base, "genes.tsv"),
                         header = FALSE, stringsAsFactors = FALSE)$V1
  rownames(mat) <- genes
  colnames(mat) <- barcodes

  # edgeR: log(CPM+1) per cell
  dge      <- DGEList(counts = mat)
  logcpm   <- cpm(dge, log = TRUE, prior.count = 1)   # log2 scale; log(CPM+1) ≈ cpm(log=T,prior=1)
  logcpm
}

message("Loading SAT counts (this may take a while)...")
sat_logcpm <- load_counts_logcpm("SAT")
message("Loading VAT counts (this may take a while)...")
vat_logcpm <- load_counts_logcpm("VAT")

# ---------- 5. load UMAP coords ----------
message("Loading UMAP coordinates...")

load_umap <- function(depot) {
  f <- file.path(base_data, depot, "ASPC_scvi_reumap_from_allobj",
                 sprintf("%s_ASPC_scVI_reUMAP.coords.tsv", depot))
  read.table(f, sep = "\t", header = TRUE, row.names = 1) %>%
    rownames_to_column("barcode") %>%
    rename(UMAP1 = 2, UMAP2 = 3)
}

sat_umap <- load_umap("SAT")
vat_umap <- load_umap("VAT")

# ---------- 6. mean log(CPM+1) for each target gene set ----------
mean_expr_on_umap <- function(logcpm_mat, umap_df, genes, set_name) {
  avail <- intersect(genes, rownames(logcpm_mat))
  message(sprintf("  %s: %d / %d genes found in matrix", set_name, length(avail), length(genes)))
  if (length(avail) == 0) return(umap_df %>% mutate(mean_expr = NA_real_))

  sub <- logcpm_mat[avail, , drop = FALSE]
  mean_vec <- if (nrow(sub) == 1) as.numeric(sub) else colMeans(sub, na.rm = TRUE)
  mean_df  <- tibble(barcode = colnames(logcpm_mat), mean_expr = mean_vec)

  umap_df %>% left_join(mean_df, by = "barcode")
}

# scale to 0-1 helper
scale01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (diff(rng) == 0) return(rep(0, length(x)))
  (x - rng[1]) / diff(rng)
}

# ---------- 7. UMAP expression plots ----------
message("Drawing UMAP expression plots...")

plot_umap_expr <- function(umap_df, title_str, scale_expr = TRUE) {
  df <- umap_df
  if (scale_expr) df <- df %>% mutate(mean_expr = scale01(mean_expr))
  expr_lim <- range(df$mean_expr, na.rm = TRUE)

  # plot cells with NA last (draw expressed cells on top)
  df <- df %>% arrange(mean_expr)

  ggplot(df, aes(x = UMAP1, y = UMAP2, color = mean_expr)) +
    geom_point(size = 0.3, alpha = 0.7) +
    scale_color_gradientn(
      colours  = expr_pal,
      limits   = expr_lim,
      na.value = "grey90",
      name     = if (scale_expr) "mean log(CPM+1)\n[0-1 scaled]" else "mean log(CPM+1)"
    ) +
    labs(title = title_str, x = "UMAP1", y = "UMAP2") +
    theme_classic(base_size = 11) +
    theme(
      plot.title     = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right",
      axis.text      = element_blank(),
      axis.ticks     = element_blank(),
      panel.border   = element_rect(colour = "black", fill = NA, linewidth = 0.8)
    ) +
    guides(color = guide_colorbar(barwidth = 0.8, barheight = 6))
}

for (kw in target_keywords) {
  genes_kw <- target_genes_combined[[kw]]
  if (length(genes_kw) == 0) {
    message(sprintf("No genes found for keyword: %s — skipping UMAP", kw))
    next
  }

  message(sprintf("Processing keyword: %s (%d genes)", kw, length(genes_kw)))

  sat_expr <- mean_expr_on_umap(sat_logcpm, sat_umap, genes_kw, paste0("SAT-", kw))
  vat_expr <- mean_expr_on_umap(vat_logcpm, vat_umap, genes_kw, paste0("VAT-", kw))

  p_sat_u <- plot_umap_expr(sat_expr, sprintf("SAT - %s (PC%d)", kw, sat_pc))
  p_vat_u <- plot_umap_expr(vat_expr, sprintf("VAT - %s (PC%d)", kw, vat_pc))

  p_combined <- p_sat_u | p_vat_u
  out_f <- file.path(out_dir,
    sprintf("SAT_PC%d.VAT_PC%d.UMAP_%s_meanExpr.pdf", sat_pc, vat_pc, kw))
  ggsave(out_f, p_combined, width = 10, height = 4.5, useDingbats = FALSE)
  message(sprintf("Saved: %s", out_f))

  # also save gene list
  write.table(
    data.frame(gene = genes_kw),
    file.path(out_dir, sprintf("%s_pathway_genes.txt", kw)),
    row.names = FALSE, quote = FALSE
  )
}

message("All done. Outputs in: ", out_dir)
