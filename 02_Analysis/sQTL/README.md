# Splicing / sQTL analysis

Splicing QTL (sQTL) analyses. Subfolders are numbered in the order they should be run:

| Folder | Description |
| --- | --- |
| `00_demultiplex/` | Prepare sample metadata, subset BAMs per donor and convert to FASTQ |
| `01_alignment/` | Alignment with STARsolo |
| `02_junction_extraction/` | Index BAMs and extract splice junctions |
| `03_meta_cell/` | Summarize cells per donor and aggregate cells into meta-cells (kNN) |
| `04_phenotype/` | Generate splicing phenotypes |
| `05_sQTL_mapping/` | sQTL mapping |
| `06_replication/` | Replication of sQTLs |
