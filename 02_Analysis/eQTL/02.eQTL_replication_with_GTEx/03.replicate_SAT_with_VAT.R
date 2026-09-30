library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
## =========================================================
## input / output
## =========================================================
data_file <- "/data/projects/11003054/e1101943/project/01.adipose_sqtl/data/09.eQTL/03.Saige_QTL/02.output/03.formal_calling/ancestry_specific_version/All_eQTL_SAT_and_VAT_with_beta.tsv"
output_dir <- "/data/projects/11003054/e1101943/project/01.adipose_sqtl/data/09.eQTL/03.Saige_QTL/02.output/04.replication/replicate_SAT_with_VAT"
cis_root <- "/data/projects/11003054/e1101943/project/01.adipose_sqtl/data/09.eQTL/03.Saige_QTL/02.output/03.formal_calling/ancestry_specific_version"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
## =========================================================
## colors
## =========================================================
colors <- c(
  "Adipocyte" = "#56BF95",
  "Early_preadipocyte"= "#c7e9c0",
  "Committed_preadipocyte" = "#a1d99b",
  "Areg" = "#CCCC4D",
  "IGFBP2+_cell" = "#74c476",
  "Mesothelial" ="#DAE220",
  "AEC" ="#F0B6DA",
  "CEC" = "#DE77AE",
  "VEC" = "#C51B7D",
  "LEC" = "#8E0052",
  "Pericyte"="#fcc5c0",
  "SMC2"="#fa9fb5",
  "SMC"="#f768a1",
  "SMC1"="#f768a1",
  "Endometrial"="#883186",
  "cDC1" = "#a6bddb",
  "cDC2" = "#74a9cf",
  "Activate_DC" = "#3690c0",
  "pDC"="#0570b0",
  "Classical_monocyte" = "#243494",
  "Non_classical_monocyte"="#061C58",
  "Intermediate_monocyte" = "#1C1894",
  "TIM4+_ATM"="#bcbddc",
  "PVM"="#9e9ac8",
  "LAM" = "#807dba",
  "TIM4+_CD11c+_ATM" = "#6a51a3",
  "B" = "#fdae61",
  "Plasmablast"="#993404",
  "Plasma_cell"="#f46d43",
  "Treg" = "#A7D5F7",
  "CD4+_T"="#9BB9FB",
  "CD8+_T"="#5087F3",
  "CD56dim_CD16+_NK" = "#4157E4",
  "CD56+_NK" = "#253CD1",
  "Mast"="#44B3E1"
)

## =========================================================
## functions
## =========================================================
get_beta_from_cis <- function(cis_root, depot, celltype, gene, marker_id) {
  cis_file <- file.path(cis_root, depot, celltype, paste0(gene, "_cis"))

  if (!file.exists(cis_file)) {
    return(NA_real_)
  }

  cis_dt <- fread(cis_file, sep = "\t", header = TRUE, data.table = TRUE, showProgress = FALSE)

  if (!("MarkerID" %in% names(cis_dt))) return(NA_real_)
  if (!("BETA" %in% names(cis_dt))) return(NA_real_)

  hit <- cis_dt[MarkerID == marker_id]

  if (nrow(hit) == 0) {
    return(NA_real_)
  } else {
    return(hit$BETA[1])
  }
}

fmt_p <- function(p) {
  ifelse(
    is.na(p), "NA",
    ifelse(p < 2.2e-16, "<2.2e-16", formatC(p, format = "e", digits = 2))
  )
}

## =========================================================
## read data
## =========================================================
data <- fread(data_file, sep = "\t", header = TRUE, data.table = TRUE)
data$gene <- gsub("-","_",data$gene)
## check required columns
required_cols <- c("depot", "celltype", "gene", "MarkerID")
missing_cols <- setdiff(required_cols, names(data))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

## if eqtl_key not present, create it
if (!("eqtl_key" %in% names(data))) {
  data[, eqtl_key := paste(gene, MarkerID, sep = "-")]
}

## only SAT and VAT
data <- data[depot %in% c("SAT", "VAT")]

## keep celltypes in a stable order
celltypes <- sort(unique(data$celltype))

all_pairs_list <- list()
stat_list <- list()

## =========================================================
## main loop
## =========================================================
stat_list <- list()

celltypes <- read.table("/data/projects/11003054/e1101943/project/01.adipose_sqtl/data/21.differential_expressed_genes/02.NEBULA/02.compare_diff_depot/Both_depot_celltypes.txt")
celltypes <- celltypes$V1

for (ct in celltypes) {
  message("Processing: ", ct)

  dt_ct <- data[celltype == ct]
  if (nrow(dt_ct) == 0) next

  sat_eqtl <- unique(dt_ct[depot == "SAT", .(gene, MarkerID)])
  vat_eqtl <- unique(dt_ct[depot == "VAT", .(gene, MarkerID)])

  sat_eqtl[, key := paste(gene, MarkerID)]
  vat_eqtl[, key := paste(gene, MarkerID)]

  sat_keys <- sat_eqtl$key
  vat_keys <- vat_eqtl$key

  ## union
  eqtl_union <- unique(dt_ct[, .(celltype, gene, MarkerID, eqtl_key)])
  if (nrow(eqtl_union) == 0) next

  ## =========================================================
  ## fetch beta
  ## =========================================================
  res_ct <- rbindlist(
    lapply(seq_len(nrow(eqtl_union)), function(i) {
      gene_i <- eqtl_union$gene[i]
      snp_i  <- eqtl_union$MarkerID[i]
      key_i  <- eqtl_union$eqtl_key[i]
      cell_i <- eqtl_union$celltype[i]

      beta_sat <- get_beta_from_cis(cis_root, "SAT", cell_i, gene_i, snp_i)
      beta_vat <- get_beta_from_cis(cis_root, "VAT", cell_i, gene_i, snp_i)

      data.table(
        celltype = cell_i,
        gene = gene_i,
        MarkerID = snp_i,
        eqtl_key = key_i,
        BETA_SAT = beta_sat,
        BETA_VAT = beta_vat
      )
    }),
    fill = TRUE
  )

  res_ct[, key := paste(gene, MarkerID)]

  res_ct[, category := fifelse(
    key %in% sat_keys & key %in% vat_keys, "Shared",
    fifelse(key %in% sat_keys & !(key %in% vat_keys), "SAT_specific",
            fifelse(!(key %in% sat_keys) & key %in% vat_keys, "VAT_specific", NA_character_))
  )]

  ## =========================================================
  ## correlation
  ## =========================================================
  stat_ct <- res_ct[, {
    dt_sub <- .SD[!is.na(BETA_SAT) & !is.na(BETA_VAT)]

    if (nrow(dt_sub) >= 3) {
      cor_test <- cor.test(dt_sub$BETA_SAT, dt_sub$BETA_VAT)
      list(
        pearson_r = unname(cor_test$estimate),
        pearson_p = cor_test$p.value,
        n = nrow(dt_sub)
      )
    } else {
      list(
        pearson_r = NA_real_,
        pearson_p = NA_real_,
        n = nrow(dt_sub)
      )
    }

  }, by = category]

  stat_ct[, celltype := ct]

  stat_list[[ct]] <- stat_ct
}

stats_dt <- rbindlist(stat_list, fill = TRUE)


fwrite(
  stats_dt,
  file = file.path(output_dir, "All_celltypes.SAT_VAT_category_correlation.tsv"),
  sep = "\t"
)
stats_dt[, category := factor(category,
                              levels = c("Shared", "SAT_specific", "VAT_specific"))]

celltype_keep <- stats_dt[, .(total_n = sum(n, na.rm = TRUE)), by = celltype][
  total_n >= 20, celltype
]

# filter
stats_dt <- stats_dt[celltype %in% celltype_keep]
p <- ggplot(stats_dt, aes(x = celltype, y = pearson_r, fill = category)) +
  geom_bar(stat = "identity",
           position = position_dodge(width = 0.8),
           width = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  scale_fill_manual(values = c(
    "Shared" = "black",
    "SAT_specific" = "grey70",
    "VAT_specific" = "grey70"
  )) +
  labs(
    x = "Cell type",
    y = "Pearson correlation (SAT vs VAT beta)",
    fill = "eQTL category",
    title = "SAT vs VAT eQTL effect size correlation by category"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(
  filename = file.path(output_dir, "SAT_VAT_correlation_barplot.pdf"),
  plot = p,
  width = 14,
  height = 6,
  useDingbats = FALSE
)

