#!/usr/bin/env Rscript
# Purpose: Compare paired donor-level immune-cell fractions between VAT and SAT.
# Input: SECTION4_METADATA, containing donor/depot and cell annotations.
# Output: SECTION4_OUTPUT/vat_sat_immune_fraction_analysis.rds.
# Method: two-sided paired t-tests with Benjamini-Hochberg correction.
# Usage: Rscript compare_vat_sat_immune_fractions.R
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
})

meta_path <- Sys.getenv("SECTION4_METADATA")
out_rds <- file.path(Sys.getenv("SECTION4_OUTPUT"), "vat_sat_immune_fraction_analysis.rds")
dir.create(dirname(out_rds), recursive = TRUE, showWarnings = FALSE)

prepare_vat_sat_ttest_data <- function(meta,
                                       cellgroup,
                                       remove_caucasian = FALSE,
                                       obesity_only = FALSE,
                                       depot_levels = c("SAT","VAT"),
                                       min_cells_per_donor_depot_ct = 10,
                                       min_total_cells_per_donor = 100,
                                       min_paired_donors_per_ct = 10) {

  squeeze01 <- function(p, eps = 1e-5) pmin(pmax(p, eps), 1 - eps)

  meta_sub <- meta %>%
    filter(Cell_group == cellgroup) %>%
    filter(!is.na(donor_id), !is.na(Depot), !is.na(Annotation_level2))

  if (remove_caucasian && "Race" %in% colnames(meta_sub)) {
    meta_sub <- meta_sub %>% filter(Race != "Caucasian")
  }

  if (obesity_only && "status" %in% colnames(meta_sub)) {
    meta_sub <- meta_sub %>% filter(status == "Obesity")
  }

  all_cts <- unique(as.character(meta_sub$Annotation_level2))
  message("\n[", cellgroup, "] Starting cell types (n=", length(all_cts), "): ",
          paste(all_cts, collapse = ", "))

  counts <- meta_sub %>%
    count(donor_id, Depot, Annotation_level2, name = "n_cells") %>%
    mutate(
      Depot = factor(Depot, levels = depot_levels),
      donor_id = factor(donor_id),
      Annotation_level2 = factor(Annotation_level2)
    )

  totals <- counts %>%
    group_by(donor_id, Depot) %>%
    summarise(total_cells = sum(n_cells), .groups = "drop")

  valid_donors <- totals %>%
    group_by(donor_id, Depot) %>%
    summarise(total_cells_all = sum(total_cells), .groups = "drop") %>%
    filter(total_cells_all >= min_total_cells_per_donor) %>%
    pull(donor_id)

  df_prop <- counts %>%
    left_join(totals, by = c("donor_id","Depot")) %>%
    mutate(
      proportion = ifelse(total_cells > 0, n_cells / total_cells, NA_real_),
      proportion = squeeze01(proportion)
    ) %>%
    filter(donor_id %in% valid_donors)

  # ── Step 1: min_cells_per_donor_depot_ct ──────────────────
  cts_before <- unique(as.character(df_prop$Annotation_level2))
  df_prop <- df_prop %>% filter(n_cells >= min_cells_per_donor_depot_ct)
  cts_after <- unique(as.character(df_prop$Annotation_level2))
  dropped <- setdiff(cts_before, cts_after)
  if (length(dropped) > 0)
    message("  [Step1: min_cells>=", min_cells_per_donor_depot_ct, "] Dropped: ", paste(dropped, collapse = ", "))

  # ── Step 2: paired depots ──────────────────────────────────
  cts_before <- unique(as.character(df_prop$Annotation_level2))
  df_prop <- df_prop %>%
    group_by(donor_id, Annotation_level2) %>%
    filter(all(depot_levels %in% Depot)) %>%
    ungroup()
  cts_after <- unique(as.character(df_prop$Annotation_level2))
  dropped <- setdiff(cts_before, cts_after)
  if (length(dropped) > 0)
    message("  [Step2: paired depots] Dropped: ", paste(dropped, collapse = ", "))

  # ── Step 3: min_paired_donors_per_ct ──────────────────────
  cts_before <- unique(as.character(df_prop$Annotation_level2))
  df_prop <- df_prop %>%
    group_by(Annotation_level2) %>%
    filter(n_distinct(donor_id) >= min_paired_donors_per_ct) %>%
    ungroup()
  cts_after <- unique(as.character(df_prop$Annotation_level2))
  dropped <- setdiff(cts_before, cts_after)
  if (length(dropped) > 0)
    message("  [Step3: min_paired_donors>=", min_paired_donors_per_ct, "] Dropped: ", paste(dropped, collapse = ", "))

  message("  [", cellgroup, "] Remaining cell types (n=", length(cts_after), "): ",
          paste(cts_after, collapse = ", "))

  if (nrow(df_prop) == 0) {
    return(list(
      ok = FALSE,
      msg = paste0("No data left after filtering. cellgroup=", cellgroup),
      cellgroup = cellgroup
    ))
  }

  ttest_results <- list()
  for (ct in unique(df_prop$Annotation_level2)) {

    df_ctw <- df_prop %>%
      filter(Annotation_level2 == ct) %>%
      select(donor_id, Depot, proportion) %>%
      distinct() %>%
      pivot_wider(names_from = Depot, values_from = proportion) %>%
      drop_na(SAT, VAT)

    if (nrow(df_ctw) < min_paired_donors_per_ct) {
      message("  [t-test skip] ", ct, ": only ", nrow(df_ctw), " paired donors after pivot")
      next
    }

    test_res <- t.test(df_ctw$VAT, df_ctw$SAT, paired = TRUE)

    ttest_results[[as.character(ct)]] <- tibble(
      Annotation_level2 = as.character(ct),
      n_paired = nrow(df_ctw),
      mean_SAT = mean(df_ctw$SAT),
      mean_VAT = mean(df_ctw$VAT),
      mean_diff = mean(df_ctw$VAT - df_ctw$SAT),
      conf.low = test_res$conf.int[1],
      conf.high = test_res$conf.int[2],
      t = unname(test_res$statistic),
      p_value = test_res$p.value
    )
  }

  df_ttest <- bind_rows(ttest_results) %>%
    mutate(
      p_adj = p.adjust(p_value, method = "fdr"),
      label = paste0("p=", signif(p_value,2), "\nFDR=", signif(p_adj,2), "\nN=", n_paired)
    )

  df_line <- df_prop %>%
    select(donor_id, Annotation_level2, Depot, proportion) %>%
    distinct() %>%
    pivot_wider(names_from = Depot, values_from = proportion) %>%
    drop_na(SAT, VAT)

  df_anno <- df_prop %>%
    group_by(Annotation_level2) %>%
    summarise(y = max(proportion, na.rm = TRUE) * 1.05, .groups = "drop") %>%
    mutate(Annotation_level2 = as.character(Annotation_level2)) %>%
    left_join(df_ttest, by = "Annotation_level2")

  list(
    ok = TRUE,
    cellgroup = cellgroup,
    settings = list(
      remove_caucasian = remove_caucasian,
      obesity_only = obesity_only,
      depot_levels = depot_levels,
      min_cells_per_donor_depot_ct = min_cells_per_donor_depot_ct,
      min_total_cells_per_donor = min_total_cells_per_donor,
      min_paired_donors_per_ct = min_paired_donors_per_ct
    ),
    df_prop = df_prop,
    df_line = df_line,
    df_ttest = df_ttest,
    df_anno = df_anno
  )
}

meta <- data.table::fread(
  meta_path, sep = "\t", data.table = FALSE,
  select = c("donor_id", "Depot", "Annotation_level1", "Annotation_level2", "Race", "status")
)

immune <- c("Mast","B","Dendritic","Macrophage","Monocyte","Neutrophil","NK","T","Plasma_cell")
Vasculars <- c("Endothelial","Pericyte","SMC")
Adipocytes <- "Adipocyte"
FAPs <- c("ASPC","Mesothelial","IGFBP2+_cell")

meta$Cell_group <- NA
meta$Cell_group[meta$Annotation_level1 %in% immune] <- "Immunes"
meta$Cell_group[meta$Annotation_level1 %in% Vasculars] <- "Vasculars"
meta$Cell_group[meta$Annotation_level1 %in% Adipocytes] <- "Adipocytes"
meta$Cell_group[meta$Annotation_level1 %in% FAPs] <- "FAPs"

data_list <- list(
  meta_path = meta_path,
  out_rds = out_rds,
  generated_at = as.character(Sys.time()),
  groups = list(
    Vasculars  = prepare_vat_sat_ttest_data(meta, "Vasculars",  remove_caucasian = TRUE, obesity_only = TRUE),
    FAPs       = prepare_vat_sat_ttest_data(meta, "FAPs",       remove_caucasian = TRUE, obesity_only = TRUE),
    Immunes    = prepare_vat_sat_ttest_data(meta, "Immunes",    remove_caucasian = TRUE, obesity_only = TRUE),
    Adipocytes = prepare_vat_sat_ttest_data(meta, "Adipocytes", remove_caucasian = TRUE, obesity_only = TRUE)
  )
)

saveRDS(data_list, out_rds)
out_rds
