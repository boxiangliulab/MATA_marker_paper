# MATA marker paper

Analysis code accompanying the MATA manuscript: a single-nucleus RNA-seq (snRNA-seq) and genetics study of human adipose tissue, covering subcutaneous (SAT) and visceral (VAT) depots across donors with different BMI, ancestry and metabolic status.

The repository contains the scripts used for data quality control, cell-type–level analyses, and genetic analyses (eQTL, sQTL, GWAS and colocalization). Raw and controlled-access data are **not** included.

## Repository structure

```
.
├── 01_Quality_control/
│   ├── Genotype/                 # Genotype QC and imputation
│   └── snRNA-seq/                # snRNA-seq processing and QC
│       ├── 01-02                 # FASTQ read modification, Cell Ranger
│       ├── 03.Demultiplexing/    # Exon variant selection, pile-up, donor demultiplexing
│       ├── 04.Library_QC/        # Library-level QC (SCTransform, demuxlet results)
│       └── 05.Integration_scanpy.py
└── 02_Analysis/
    ├── CellProportionAnalysis/             # Cell-type proportion analysis
    ├── Differential_expression_gene_analysis/  # Pseudobulk DE, BMI x depot interaction (limma)
    ├── Adipogenesis_in_committedASPC/      # Committed ASPC proportions, module scores, GSEA
    ├── SMC_phenotypic_change/              # Smooth muscle cell contractile/synthetic switching
    ├── VAT_DP_ATM_analysis/                # VAT DP adipose tissue macrophage analyses
    ├── eQTL/                               # eQTL mapping, GTEx replication, snpEff, TSS,
    │                                       # cell-state & dynamic eQTL, heritability, coloc
    ├── sQTL/                               # Splicing QTL: alignment, junctions, meta-cells,
    │                                       # phenotypes, mapping, replication
    ├── GWAS/                               # GWAS harmonization, liftover to hg38, allele alignment
    └── Colocalization/                     # Colocalization between xQTLs and GWAS
```

## Workflow overview

1. **Genotype QC** (`01_Quality_control/Genotype`) – pre-imputation QC, imputation, and post-imputation QC.
2. **snRNA-seq QC** (`01_Quality_control/snRNA-seq`) – read processing, Cell Ranger alignment, genotype-based donor demultiplexing, library QC and integration with Scanpy.
3. **Cell-level analyses** (`02_Analysis`) – cell-type proportions, pseudobulk differential expression, adipogenesis in committed ASPCs, SMC phenotypic change, and VAT DP-ATM analyses.
4. **Genetic analyses** (`02_Analysis/eQTL`, `sQTL`, `GWAS`, `Colocalization`) – QTL mapping and replication, GWAS harmonization, and colocalization.

Scripts within each folder are numbered in the order they should be run. Some subfolders contain their own `README.md` with more detail (e.g. `02_Analysis/VAT_DP_ATM_analysis`).

## Requirements

The code uses **R**, **Python** and **shell** scripts. Main tools include Cell Ranger, demuxlet, Scanpy, scVelo, CellPhoneDB, Seurat/SCTransform, limma and msigdbr. A conda environment file is provided for the VAT DP-ATM analyses (`02_Analysis/VAT_DP_ATM_analysis/environment.yml`).

## Data availability

Input data (raw sequencing data, genotypes, and cell-level metadata) are not distributed with this repository. Paths in the scripts need to be adapted to your local environment.

## Citation

If you use this code, please cite the MATA manuscript (citation to be added upon publication).

## Contact

Liu Boxiang Lab — please open an issue in this repository for questions.
