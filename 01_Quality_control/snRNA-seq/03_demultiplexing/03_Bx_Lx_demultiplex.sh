#PBS -q normal
#PBS -l select=1:ncpus=32:mem=200G
#PBS -l walltime=2:00:00
#PBS -N Bx_Lx_demultiplex
#PBS -o Bx_Lx_demultiplex.o
#PBS -e Bx_Lx_demultiplex.e
#PBS -P 11003054
#PBS -j oe

cd /project/adipose/demuxlet/02.demultiplex

if [ ! -d "Groupx" ]; then
  mkdir Groupx
fi
cd Groupx
mkdir -p Bx_Lx
cd Bx_Lx

.conda/envs/sr/bin/popscle demuxlet \
  --plp /project/adipose/demuxlet/01.pileup/Groupx/Bx_Lx/Bx_Lx \
  --group-list /project/adipose/demuxlet/00.renamed_barcode/barcodes_Bx_Lx.tsv.gz \
  --doublet-prior 0.25 \
  --geno-error-coeff 0.1 \
  --field GT \
  --vcf /project/adipose/genotype/02.vcf_bygroup/Groupx_GT.vcf \
  --out Bx_Lx
