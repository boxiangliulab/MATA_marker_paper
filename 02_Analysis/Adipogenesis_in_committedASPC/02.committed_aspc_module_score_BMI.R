# Module scores for committed preadipocytes ~ BMI
#
# Gene sets (from 06c/06d):
#   Collagen/fibrosis:
#     Eisinger2023, NABA_COLLAGENS, REACTOME_COLLAGEN_FORMATION
#   Adipogenesis:
#     TF_adipogenesis, HALLMARK_ADIPOGENESIS, REACTOME_ADIPOGENESIS,
#     GOBP_WHITE_FAT_CELL_DIFFERENTIATION, KEGG_PPAR_SIGNALING_PATHWAY
#
# Method:
#   Seurat AddModuleScore (per-cell) -> mean per donor per depot
#   LM: score ~ BMI_group_num + Age + Gender + T2DM + Race
#   Spearman: score vs BMI_group_num
#
# Input:
#   committed_aspc_counts/{SAT,VAT}/  (10x MEX from 00.export_committed_aspc_counts.py)
# Output:
#   02.committed_aspc_module_score_BMI_boxplot.pdf
#   02.committed_aspc_module_score_BMI_loess.pdf
#   02.committed_aspc_module_score_BMI_stats.tsv

suppressPackageStartupMessages({
  library(Seurat)
  library(msigdbr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(data.table)
})

IN_DIR    <- "/data/projects/11003054/e1101943/project/01.adipose_sqtl/data/30.adipogenesis/committed_aspc_counts"
META_FILE <- "/data/projects/11003054/e1101943/project/01.adipose_sqtl/data/08.integration_snRNA_object/01.freeze_version/metadata_after_libraryQC_after_ctQC_freeze.tsv"
OUT_DIR   <- "/data/projects/11003054/e1101943/project/01.adipose_sqtl/data/30.adipogenesis"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

COL_VAT <- "#F3AF44"
COL_SAT <- "#67B3AD"
DEPOT_COLORS <- c(SAT = COL_SAT, VAT = COL_VAT)

GROUP_COLORS <- c(
  "Non-obese"      = "#efedf5",
  "Obese"          = "#bcbddc",
  "Severely obese" = "#807dba",
  "Morbidly obese" = "#54278f"
)

ASIAN_RACES <- c("Chinese", "Malay", "Indian")

# ── Gene sets ─────────────────────────────────────────────────────────────────
cat("Fetching MSigDB gene sets...\n")
msig_all <- msigdbr(species = "Homo sapiens")
get_msig <- function(gs) unique(msig_all$gene_symbol[msig_all$gs_name == gs])

GENESETS <- list(
  # Collagen / fibrosis (06c)
  Eisinger2023 = c(
    "COL1A1", "COL3A1", "COL6A1", "COL6A3", "LUM", "TNC",
    "COL5A1", "LOX", "FN1", "LOXL1", "LOXL2", "CCN2",
    "TGFB1", "VCAN", "LRP1", "THBS1", "THY1", "ACTA2"
  ),
  NABA_Collagens            = get_msig("NABA_COLLAGENS"),
  Collagen_Formation        = get_msig("REACTOME_COLLAGEN_FORMATION"),
  # Adipogenesis (06d)
  HALLMARK_ADIPOGENESIS     = get_msig("HALLMARK_ADIPOGENESIS"),
  REACTOME_ADIPOGENESIS     = get_msig("REACTOME_ADIPOGENESIS"),
  GOBP_WFC                  = get_msig("GOBP_WHITE_FAT_CELL_DIFFERENTIATION"),
  KEGG_PPAR                 = get_msig("KEGG_PPAR_SIGNALING_PATHWAY"))

for (nm in names(GENESETS))
  cat(sprintf("  %-35s %d genes\n", nm, length(GENESETS[[nm]])))

# ── 1. Metadata ───────────────────────────────────────────────────────────────
meta_full <- fread(META_FILE, data.table = FALSE, na.strings = c("", "NA"))

donor_meta <- meta_full %>%
  filter(Race %in% ASIAN_RACES, !is.na(BMI)) %>%
  mutate(
    BMI           = as.numeric(BMI),
    T2DM          = ifelse(Diabetic == "T2DM", 1, 0),
    BMI_group     = case_when(
      BMI <  27.5 ~ "Non-obese",
      BMI <  32.5 ~ "Obese",
      BMI <  37.5 ~ "Severely obese",
      TRUE        ~ "Morbidly obese"
    ),
    BMI_group     = factor(BMI_group,
                           levels = c("Non-obese", "Obese",
                                      "Severely obese", "Morbidly obese")),
    BMI_group_num = as.numeric(BMI_group)
  ) %>%
  select(donor_id, BMI, BMI_group, BMI_group_num, Age, Gender, Race, T2DM) %>%
  distinct(donor_id, .keep_all = TRUE)

# ── 2. Load 10x -> Seurat, normalize, add module scores ──────────────────────
build_seurat <- function(depot) {
  cat(sprintf("Loading %s ...\n", depot))
  mat <- Read10X(file.path(IN_DIR, depot))
  seu <- CreateSeuratObject(counts = mat, project = depot, min.cells = 0, min.features = 0)
  seu <- NormalizeData(seu, verbose = FALSE)

  # keep only cells whose donor is in donor_meta
  # cell barcodes contain donor_id (format: barcode-donor_id or similar)
  # safer: pull donor_id from exported metadata
  meta_path <- file.path(IN_DIR, depot, "metadata.tsv.gz")
  if (file.exists(meta_path)) {
    cell_meta <- read.table(meta_path, sep = "\t", header = TRUE,
                            row.names = 1, check.names = FALSE)
    # align to Seurat cell names
    common <- intersect(rownames(cell_meta), colnames(seu))
    seu     <- seu[, common]
    cell_meta <- cell_meta[common, , drop = FALSE]
    seu <- AddMetaData(seu, cell_meta)
  }

  seu$Depot <- depot
  seu
}

seu_sat <- build_seurat("SAT")
seu_vat <- build_seurat("VAT")

# AddModuleScore per depot (avoids gene universe mismatch)
run_module_scores <- function(seu) {
  gs_list <- lapply(GENESETS, function(g) intersect(g, rownames(seu)))
  # remove sets with < 3 genes
  gs_list <- gs_list[sapply(gs_list, length) >= 3]
  cat(sprintf("  [%s] computing %d module scores\n", unique(seu$Depot), length(gs_list)))
  seu <- AddModuleScore(seu, features = gs_list, name = "MS_", seed = 42)
  # rename columns: MS_1 -> gs_name
  score_cols <- paste0("MS_", seq_along(gs_list))
  for (i in seq_along(gs_list))
    seu@meta.data[[names(gs_list)[i]]] <- seu@meta.data[[score_cols[i]]]
  seu@meta.data[score_cols] <- NULL
  seu
}

seu_sat <- run_module_scores(seu_sat)
seu_vat <- run_module_scores(seu_vat)

# ── 3. Aggregate per donor (mean module score per donor per depot) ────────────
aggregate_scores <- function(seu, depot_name) {
  md <- seu@meta.data
  gs_cols <- names(GENESETS)[names(GENESETS) %in% colnames(md)]
  if (!"donor_id" %in% colnames(md))
    stop("donor_id not found in Seurat metadata for ", depot_name)
  md %>%
    group_by(donor_id) %>%
    summarise(across(all_of(gs_cols), mean, na.rm = TRUE), .groups = "drop") %>%
    mutate(Depot = depot_name)
}

scores_sat <- aggregate_scores(seu_sat, "SAT")
scores_vat <- aggregate_scores(seu_vat, "VAT")

score_df <- bind_rows(scores_sat, scores_vat) %>%
  left_join(donor_meta, by = "donor_id") %>%
  filter(!is.na(BMI)) %>%
  mutate(Depot = factor(Depot, levels = c("SAT", "VAT")))

# long format
gs_cols  <- names(GENESETS)[names(GENESETS) %in% colnames(score_df)]
score_long <- score_df %>%
  pivot_longer(all_of(gs_cols), names_to = "gs_name", values_to = "score") %>%
  mutate(gs_name = factor(gs_name, levels = names(GENESETS)))

# ── 3b. Save per-cell module scores ──────────────────────────────────────────
cell_scores <- bind_rows(
  seu_sat@meta.data %>% mutate(cell_barcode = rownames(.)),
  seu_vat@meta.data %>% mutate(cell_barcode = rownames(.))
) %>%
  select(cell_barcode, donor_id, Depot, any_of(names(GENESETS)))

fwrite(cell_scores,
       file.path(OUT_DIR, "02.committed_aspc_module_score_per_cell.tsv.gz"),
       sep = "\t", quote = FALSE, row.names = FALSE)
cat("Saved: 02.committed_aspc_module_score_per_cell.tsv.gz\n")

# ── 4. LM + Spearman ─────────────────────────────────────────────────────────
fmt_sci     <- function(p) formatC(p, format = "e", digits = 2)
safe_covars <- function(df) {
  covars <- c("Age", "Gender", "T2DM", "Race")
  covars[sapply(covars, function(v) {
    x <- df[[v]]
    if (is.numeric(x)) sd(x, na.rm = TRUE) > 0
    else length(unique(na.omit(x))) >= 2
  })]
}

stats_res <- score_long %>%
  group_by(gs_name, Depot) %>%
  group_modify(~ {
    keep  <- safe_covars(.x)
    fmla  <- as.formula(paste("score ~ BMI_group_num +", paste(keep, collapse = " + ")))
    m     <- lm(fmla, data = .x)
    ct    <- summary(m)$coefficients
    sp    <- cor.test(.x$score, .x$BMI_group_num, method = "spearman", exact = FALSE)
    data.frame(
      n_donors = nrow(.x),
      lm_beta  = ct["BMI_group_num", "Estimate"],
      lm_p     = ct["BMI_group_num", "Pr(>|t|)"],
      sp_rho   = sp$estimate,
      sp_p     = sp$p.value
    )
  }) %>% ungroup() %>%
  mutate(lm_p_adj = p.adjust(lm_p, method = "BH"),
         sp_p_adj = p.adjust(sp_p, method = "BH"))

write.table(stats_res,
            file.path(OUT_DIR, "02.committed_aspc_module_score_BMI_stats.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("\n=== LM results ===\n")
print(stats_res %>% select(gs_name, Depot, n_donors, lm_beta, lm_p, lm_p_adj, sp_rho))

# ── 4b. Interaction term: BMI_group_num × Depot ───────────────────────────────
# Combined SAT+VAT model to test whether BMI effect differs between depots
safe_covars_combined <- function(df) {
  covars <- c("Age", "Gender", "T2DM", "Race")
  covars[sapply(covars, function(v) {
    x <- df[[v]]
    if (is.numeric(x)) sd(x, na.rm = TRUE) > 0
    else length(unique(na.omit(x))) >= 2
  })]
}

interaction_res <- score_long %>%
  group_by(gs_name) %>%
  group_modify(~ {
    df <- .x %>% mutate(Depot = factor(Depot, levels = c("SAT", "VAT")))
    keep <- safe_covars_combined(df)
    rhs  <- paste(c("BMI_group_num * Depot", keep), collapse = " + ")
    fmla <- as.formula(paste("score ~", rhs))
    m    <- tryCatch(lm(fmla, data = df), error = function(e) NULL)
    if (is.null(m)) return(data.frame(n_SAT=NA, n_VAT=NA,
                                       int_beta=NA, int_p=NA))
    ct <- summary(m)$coefficients
    int_term <- "BMI_group_num:DepotVAT"
    if (!int_term %in% rownames(ct))
      return(data.frame(n_SAT=NA, n_VAT=NA, int_beta=NA, int_p=NA))
    data.frame(
      n_SAT    = sum(df$Depot == "SAT"),
      n_VAT    = sum(df$Depot == "VAT"),
      int_beta = ct[int_term, "Estimate"],
      int_p    = ct[int_term, "Pr(>|t|)"]
    )
  }) %>% ungroup() %>%
  mutate(
    int_p_adj_BH        = p.adjust(int_p, method = "BH"),
    int_p_adj_bonf      = p.adjust(int_p, method = "bonferroni"),
    int_p_adj_holm      = p.adjust(int_p, method = "holm"),
    int_p_adj_BY        = p.adjust(int_p, method = "BY")
  )

# q-value (Storey) — only if qvalue package available
if (requireNamespace("qvalue", quietly=TRUE)) {
  qobj <- tryCatch(qvalue::qvalue(interaction_res$int_p[!is.na(interaction_res$int_p)]),
                   error=function(e) NULL)
  if (!is.null(qobj)) {
    interaction_res$int_qvalue <- NA_real_
    interaction_res$int_qvalue[!is.na(interaction_res$int_p)] <- qobj$qvalues
  }
}

cat("\n=== Interaction: BMI_group_num × Depot ===\n")
cat("(int_beta = VAT slope minus SAT slope; positive = VAT stronger BMI effect)\n")
cat("\nRaw p-values:\n")
print(interaction_res %>% arrange(int_p) %>%
        dplyr::select(gs_name, int_beta, int_p,
                      int_p_adj_BH, int_p_adj_bonf, int_p_adj_holm, int_p_adj_BY,
                      matches("int_qvalue")))

write.table(interaction_res,
            file.path(OUT_DIR, "02.committed_aspc_module_score_interaction_stats.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ── 5. Boxplot (BMI group, facet by depot, one panel per gene set) ────────────
n_labels_box <- score_long %>%
  distinct(donor_id, Depot, BMI_group) %>%
  count(Depot, BMI_group)

box_plots <- lapply(levels(score_long$gs_name), function(gs) {
  df  <- score_long %>% filter(gs_name == gs)
  sdf <- stats_res  %>% filter(gs_name == gs) %>%
    mutate(label = paste0(
      "beta=", round(lm_beta, 3), ", p=", fmt_sci(lm_p),
      ifelse(lm_p_adj < 0.05, " *", "")
    ))

  ggplot(df, aes(x = BMI_group, y = score, fill = BMI_group)) +
    geom_boxplot(outlier.size = 0.4, width = 0.37, color = "grey30", linewidth = 0.3) +
    geom_jitter(width = 0.15, size = 0.3, color = "grey40", alpha = 0.5) +
    geom_smooth(aes(x = BMI_group_num, y = score),
                method = "lm", formula = y ~ x,
                color = "red", fill = "red", alpha = 0.15,
                linewidth = 0.5, inherit.aes = FALSE) +
    geom_text(data = n_labels_box,
              aes(x = BMI_group, y = -Inf, label = paste0("n=", n)),
              inherit.aes = FALSE, vjust = -0.3, size = 1.9, color = "grey40") +
    geom_text(data = sdf,
              aes(x = 2.5, y = Inf, label = label),
              inherit.aes = FALSE, vjust = 1.1, hjust = 0.5, size = 2.0, lineheight = 0.9) +
    facet_wrap(~ Depot, nrow = 1) +
    scale_fill_manual(values = GROUP_COLORS) +
    labs(x = NULL, y = "Module score",
         title = ifelse(!is.na(MODULE_LABELS[gs]), MODULE_LABELS[gs], gsub("_", " ", gs))) +
    theme_classic(base_size = 9) +
    theme(legend.position  = "none",
          axis.text.x      = element_text(size = 7, angle = 90, vjust = 0.5, hjust = 1),
          strip.background = element_blank(),
          strip.text       = element_text(face = "bold", size = 8),
          plot.title       = element_text(size = 8, face = "bold"))
})

ggsave(file.path(OUT_DIR, "02.committed_aspc_module_score_BMI_boxplot.pdf"),
       wrap_plots(box_plots, ncol = 2),
       width = 5.5, height = length(box_plots) * 1.3, limitsize = FALSE)
cat("Saved: 02.committed_aspc_module_score_BMI_boxplot.pdf\n")

# ── 6. Heatmap: mean z-score per module × BMI group × depot ──────────────────
# rows = module (fibrosis then adipogenesis), cols = BMI group
# facet by depot; color = mean z-score; stars = LM significance

GS_GROUPS <- list(
  Fibrosis     = c("Eisinger2023", "NABA_Collagens", "Collagen_Formation"),
  Adipogenesis = c("HALLMARK_ADIPOGENESIS", "REACTOME_ADIPOGENESIS",
                   "GOBP_WFC", "KEGG_PPAR")
)
MODULE_ORDER <- c(rev(GS_GROUPS$Fibrosis), rev(GS_GROUPS$Adipogenesis))

MODULE_LABELS <- c(
  Eisinger2023            = "18 literature\nfibrosis genes",
  NABA_Collagens          = "NABA Collagens",
  Collagen_Formation      = "Collagen Formation",
  HALLMARK_ADIPOGENESIS   = "HALLMARK Adipogenesis",
  REACTOME_ADIPOGENESIS   = "REACTOME Adipogenesis",
  GOBP_WFC                = "GOBP White Fat Cell Diff.",
  KEGG_PPAR               = "KEGG PPAR Signaling"
)

# z-score each module × depot across donors
score_z <- score_long %>%
  group_by(gs_name, Depot) %>%
  mutate(score_z = as.numeric(scale(score))) %>%
  ungroup()

# mean z-score per module × depot × BMI group
heat_df <- score_z %>%
  group_by(gs_name, Depot, BMI_group) %>%
  summarise(mean_z = mean(score_z, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    gs_name   = factor(gs_name,   levels = MODULE_ORDER),
    BMI_group = factor(BMI_group, levels = c("Non-obese", "Obese",
                                             "Severely obese", "Morbidly obese")),
    Depot     = factor(Depot, levels = c("SAT", "VAT"))
  )

# significance stars from LM
star_label <- function(p) case_when(
  p < 0.001 ~ "***",
  p < 0.01  ~ "**",
  p < 0.05  ~ "*",
  p < 0.1   ~ ".",
  TRUE       ~ ""
)

sig_df <- stats_res %>%
  mutate(
    stars   = star_label(lm_p),
    gs_name = factor(gs_name, levels = MODULE_ORDER),
    Depot   = factor(Depot, levels = c("SAT", "VAT"))
  ) %>%
  # place star at last BMI group column
  mutate(BMI_group = factor("Morbidly obese",
                            levels = levels(heat_df$BMI_group)))

# category strip (fibrosis vs adipogenesis separator)
category_df <- tibble(
  gs_name  = factor(MODULE_ORDER, levels = MODULE_ORDER),
  category = c(rep("Fibrosis", 3), rep("Adipogenesis", 4))
)

HEAT_COLORS <- c("#313695","#4575b4","#74add1","#abd9e9","#e0f3f8",
                 "#ffffbf","#fee090","#fdae61","#f46d43","#d73027","#a50026")

p_heat <- ggplot(heat_df, aes(x = BMI_group, y = gs_name, fill = mean_z)) +
  geom_tile(color = "white", linewidth = 0.4) +
  # significance stars
  geom_text(data = sig_df,
            aes(x = BMI_group, y = gs_name, label = stars),
            inherit.aes = FALSE, size = 3, vjust = 0.75) +
  facet_wrap(~ Depot, nrow = 1) +
  scale_fill_gradientn(
    colors   = HEAT_COLORS,
    name     = "Mean z-score",
    na.value = "grey90"
  ) +
  scale_y_discrete(labels = MODULE_LABELS) +
  # horizontal separator between fibrosis and adipogenesis
  geom_hline(yintercept = 3.5, color = "grey50", linewidth = 0.6, linetype = "dashed") +
  annotate("text", x = 0.4, y = 5.5, label = "Adipogenesis",
           hjust = 0, size = 2.5, color = "grey30", fontface = "italic") +
  annotate("text", x = 0.4, y = 2,   label = "Fibrosis",
           hjust = 0, size = 2.5, color = "grey30", fontface = "italic") +
  labs(x = NULL, y = NULL,
       title = "Committed ASPC module scores ~ BMI group",
       caption = "Stars: LM p < 0.1 (.), < 0.05 (*), < 0.01 (**), < 0.001 (***)") +
  theme_classic(base_size = 9) +
  theme(
    axis.text.x      = element_text(size = 7, angle = 30, vjust = 1, hjust = 1),
    axis.text.y      = element_text(size = 8),
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = 9),
    plot.title       = element_text(size = 9, face = "bold"),
    plot.caption     = element_text(size = 7, color = "grey50"),
    legend.key.height = unit(0.8, "cm"),
    legend.text      = element_text(size = 7),
    legend.title     = element_text(size = 8)
  )

ggsave(file.path(OUT_DIR, "02.committed_aspc_module_score_heatmap.pdf"),
       p_heat, width = 4.5, height = 3.2)
cat("Saved: 02.committed_aspc_module_score_heatmap.pdf\n")


cat("\nDone.\n")
