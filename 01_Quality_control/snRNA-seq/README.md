# snRNA-seq QC

Processing and quality control of the snRNA-seq data. Run in order:

| Step | Description |
| --- | --- |
| `01.Modify_fastq_read.sh` | Modify FASTQ reads before alignment |
| `02.CellRanger.sh` | Alignment and counting with Cell Ranger |
| `03.Demultiplexing/` | Select exon variant sites, merge VCFs, pile-up and genotype-based donor demultiplexing (demuxlet) |
| `04.Library_QC/` | Library-level QC with SCTransform, using demuxlet results (with and without genotypes) |
| `05.Integration_scanpy.py` | Data integration with Scanpy |
