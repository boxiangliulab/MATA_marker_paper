# VAT DP ATM analysis

This directory contains analysis scripts supporting the VAT DP-ATM analyses in
Section 4 of the MATA manuscript. It intentionally excludes figure assembly,
panel formatting, microscopy processing, plotting-ready RDS construction, and
the drug-target analysis removed from the final manuscript.

## Workflow

1. `00_pseudobulk/`: aggregate cell-level raw counts into donor-level count and
   TMM-normalized expression matrices, and prepare VAT donor groups and cell
   counts used by downstream analyses.
2. `01_macrophage_subtypes/`: compare VAT and SAT immune-cell fractions,
   identify macrophage subtype markers, and perform pathway enrichment.
3. `02_rna_velocity/`: fit the SP-to-DP ATM scVelo dynamical model.
4. `03_cell_communication/`: run VAT CellPhoneDB analysis and extract
   DP-ATM IL1B/IL1-receptor interactions.
5. `04_module_scores/`: calculate donor-level inflammatory and insulin-receptor
   signaling module scores.
6. `05_mesothelial_associations/`: test DP-ATM associations with mesothelial
   fibrosis genes, cytokines, and insulin-receptor signaling.
7. `06_clinical_associations/`: test SP/DP abundance associations with BMI,
   T2D, and HbA1c.
8. `07_sensitivity_analyses/`: repeat key models with ATM abundance and capture
   depth covariates or continuous DP abundance.

## Input data

Large or controlled-access inputs are not included. Each script documents its
required inputs and outputs in its header. Expected inputs include:

- frozen cell-level metadata with donor, depot, ancestry, clinical variables,
  and level-1/level-2 annotations;
- cell-type raw-count matrices;
- annotated macrophage H5AD and velocyto loom files;
- CellPhoneDB v5 reference data;
- MSigDB gene sets obtained through `msigdbr`.

Example paths are listed in `config/paths.example.env`. No script should contain
user-specific or server-specific absolute paths.

Copy the example file, replace every placeholder, and export the variables
before running the scripts:

```bash
cp config/paths.example.env config/paths.env
set -a
source config/paths.env
set +a
```

Run scripts from the repository root in the numbered workflow order. Individual
scripts can be rerun once their listed inputs exist. `SCRIPT_INDEX.tsv` provides
a compact input/output map for every included script.

## Analysis conventions

- Analyses are restricted to VAT unless explicitly stated otherwise.
- The macrophage-subtype DEG and GSEA analysis combines SAT and VAT
  macrophages. For each subtype, Scanpy performs a one-versus-rest comparison;
  the DP result therefore compares DP ATMs with LAM, PVM, and SP ATMs combined.
- DP ATM denotes `TIM4+_CD11c+_ATM`; SP ATM denotes `TIM4+_ATM`.
- DP-present donors have at least 10 DP ATM nuclei; DP-absent donors have zero.
  Donors with 1-9 DP ATM nuclei are excluded from binary DP-status analyses.
- Donors require at least 10 nuclei for donor-level cell-type summaries unless
  a script states a different threshold.
- Caucasian donors are excluded from the Section 4 association analyses.
- Donors with zero mesothelial IL6 or CXCL8 expression are excluded from the
  corresponding cytokine association analysis.
- The SP-DP fraction analysis reports the raw regression slope and the
  covariate-adjusted p-value. The T2D analysis reports both the unadjusted
  chi-square test and the age-, gender-, race-, and BMI-adjusted odds ratio.

## Outputs

Scripts write analysis tables and intermediate analytical objects to a
user-specified output directory. Generated data and manuscript figures are not
stored in this source-code directory.
