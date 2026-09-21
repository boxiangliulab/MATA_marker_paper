#!/usr/bin/env Rscript
# Purpose: Generate the donor-level VAT group and cell-count inputs used by the
# DP-status, mesothelial, and sensitivity analyses.
# Input: SECTION4_METADATA.
# Outputs: VAT_donor_groups.tsv and VAT_ATM_counts_by_donor.tsv under
# SECTION4_OUTPUT/donor_summaries.
# Method: DP- and SP-ATM counts are calculated from level-2 annotations. Donors
# with zero target cells are absent, donors with at least 10 are present, and
# donors with 1-9 are assigned NA for the corresponding binary comparison.
# Usage: Rscript prepare_vat_donor_summaries.R

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
})

metadata_path <- Sys.getenv("SECTION4_METADATA")
output_dir <- file.path(Sys.getenv("SECTION4_OUTPUT"), "donor_summaries")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

metadata <- fread(metadata_path, sep = "\t", data.table = FALSE)
vat <- metadata %>% filter(Depot == "VAT")

donor_meta_columns <- c(
  "donor_id", "Age", "Gender", "Race", "BMI", "status", "Self_defined_DM",
  "b_FBG", "b_HbA1c", "b_Insulin (mU/L)", "Body Fat (%)"
)
missing_columns <- setdiff(
  c(donor_meta_columns, "Annotation_level1", "Annotation_level2"),
  colnames(vat)
)
if (length(missing_columns) > 0) {
  stop("Missing metadata columns: ", paste(missing_columns, collapse = ", "))
}

donor_meta <- vat %>%
  select(all_of(donor_meta_columns)) %>%
  distinct(donor_id, .keep_all = TRUE)

dp_counts <- vat %>%
  filter(Annotation_level2 == "TIM4+_CD11c+_ATM") %>%
  count(donor_id, name = "n_DP")

sp_counts <- vat %>%
  filter(Annotation_level2 == "TIM4+_ATM") %>%
  count(donor_id, name = "n_TIM")

donor_groups <- donor_meta %>%
  left_join(dp_counts, by = "donor_id") %>%
  left_join(sp_counts, by = "donor_id") %>%
  mutate(
    n_DP = replace_na(n_DP, 0L),
    n_TIM = replace_na(n_TIM, 0L),
    group_DP = case_when(
      n_DP >= 10 ~ "DP_present",
      n_DP == 0 ~ "DP_absent",
      TRUE ~ NA_character_
    ),
    group_TIM = case_when(
      n_TIM >= 10 ~ "TIM4_present",
      n_TIM == 0 ~ "TIM4_absent",
      TRUE ~ NA_character_
    )
  )

vat_cell_counts <- vat %>% count(donor_id, name = "n_VAT_cells")
atm_counts <- vat %>%
  filter(Annotation_level1 == "Macrophage") %>%
  count(donor_id, name = "n_ATM")

donor_counts <- donor_meta %>%
  select(donor_id) %>%
  left_join(atm_counts, by = "donor_id") %>%
  left_join(dp_counts, by = "donor_id") %>%
  left_join(vat_cell_counts, by = "donor_id") %>%
  mutate(across(c(n_ATM, n_DP, n_VAT_cells), ~ replace_na(.x, 0L)))

groups_path <- file.path(output_dir, "VAT_donor_groups.tsv")
counts_path <- file.path(output_dir, "VAT_ATM_counts_by_donor.tsv")
fwrite(donor_groups, groups_path, sep = "\t")
fwrite(donor_counts, counts_path, sep = "\t")

message("[saved] ", groups_path)
message("[saved] ", counts_path)
