#!/usr/bin/env Rscript
# Purpose: Test SP- and DP-ATM fractions against BMI in VAT.
# Input: SECTION4_METADATA.
# Output: SECTION4_OUTPUT/sp_dp_fraction_vs_bmi.tsv.
# Method: Fractions are relative to all VAT macrophages. For each subtype,
# donors with at least 10 cells are retained. The model is:
#   subtype fraction ~ BMI + Age + Gender + Race + Self_defined_DM
# Caucasian donors are excluded.
# Usage: Rscript sp_dp_fraction_bmi_association.R

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

metadata <- Sys.getenv("SECTION4_METADATA")
output <- file.path(Sys.getenv("SECTION4_OUTPUT"), "sp_dp_fraction_vs_bmi.tsv")
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)

meta <- fread(
  metadata,
  select = c("donor_id", "Depot", "Annotation_level1", "Annotation_level2",
             "Age", "Gender", "Race", "BMI", "Self_defined_DM"),
  data.table = FALSE
) %>%
  filter(Depot == "VAT", Race != "Caucasian")

covariates <- meta %>%
  group_by(donor_id) %>%
  summarise(
    across(c(Age, Gender, Race, BMI, Self_defined_DM), ~ first(na.omit(.x))),
    .groups = "drop"
  ) %>%
  mutate(
    Age = suppressWarnings(as.numeric(as.character(Age))),
    BMI = suppressWarnings(as.numeric(as.character(BMI))),
    Gender = factor(Gender),
    Race = factor(Race),
    Self_defined_DM = factor(if_else(Self_defined_DM == "T2DM", "T2DM", "Non_T2DM"))
  )

macrophage_counts <- meta %>%
  filter(Annotation_level1 == "Macrophage") %>%
  count(donor_id, name = "n_macrophage")

subtypes <- c(SP_ATM = "TIM4+_ATM", DP_ATM = "TIM4+_CD11c+_ATM")

results <- lapply(names(subtypes), function(label) {
  subtype_counts <- meta %>%
    filter(Annotation_level2 == subtypes[[label]]) %>%
    count(donor_id, name = "n_subtype") %>%
    filter(n_subtype >= 10)

  analysis_data <- macrophage_counts %>%
    inner_join(subtype_counts, by = "donor_id") %>%
    mutate(fraction = n_subtype / n_macrophage) %>%
    inner_join(covariates, by = "donor_id") %>%
    filter(complete.cases(.))

  fit <- lm(fraction ~ BMI + Age + Gender + Race + Self_defined_DM,
            data = analysis_data)
  stat <- coef(summary(fit))["BMI", ]

  data.frame(
    subtype = label,
    n_donors = nrow(analysis_data),
    beta_BMI = unname(stat["Estimate"]),
    standard_error = unname(stat["Std. Error"]),
    p_value = unname(stat["Pr(>|t|)"])
  )
}) %>% bind_rows()

fwrite(results, output, sep = "\t")
print(results)
message("[saved] ", output)
