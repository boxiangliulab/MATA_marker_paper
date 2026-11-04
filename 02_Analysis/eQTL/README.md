# eQTL analysis

Expression QTL (eQTL) analyses. Subfolders are numbered in the order they should be run:

| Folder | Description |
| --- | --- |
| `01.eQTL_mapping/` | Export counts from Scanpy, generate phenotypes and run eQTL mapping (SAIGE-QTL) |
| `02.eQTL_replication_with_GTEx/` | Replication with GTEx and between SAT and VAT |
| `03.snpEff/` | Variant annotation with snpEff and functional enrichment |
| `04.tss_region/` | Distance of lead SNPs to the TSS |
| `05.cell_state_eQTL/` | Cell-state eQTL analysis |
| `06.Dynamic_eQTL/` | Dynamic eQTL analysis |
| `07.Heritability/` | Heritability analysis |
| `08.colocaliztion/` | Colocalization of eQTLs with GWAS |
