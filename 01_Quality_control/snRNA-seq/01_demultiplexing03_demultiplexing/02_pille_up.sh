#PBS -q normal
#PBS -l select=1:ncpus=32:mem=400gb
#PBS -l walltime=24:00:00
#PBS -N Bx_Lx_PE
#PBS -o Bx_Lx_PE.o
#PBS -e Bx_Lx_PE.e
#PBS -P 11003054


cd /project/adipose/demuxlet/01.pileup

if [ ! -d "Groupx" ]; then
  mkdir Groupx
fi
cd Groupx
mkdir -p Bx_Lx
cd Bx_Lx

.conda/envs/sr/bin/popscle dsc-pileup \
  --sam /project/01.adipose_sqtl/data/07.full_scale_snRNA/04.cellranger/Groupx/Bx_Lx/outs/possorted_genome_bam.bam \
  --vcf /project/adipose/genotype/02.vcf_bygroup/Groupx_GT.vcf \
  --group-list /project/adipose/demuxlet/00.renamed_barcode/barcodes_Bx_Lx.tsv.gz \
  --out Bx_Lx
