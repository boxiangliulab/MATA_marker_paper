#!/bin/bash
#PBS -q normal
#PBS -l select=1:ncpus=2:mem=500G
#PBS -l walltime=2:00:00
#PBS -N 16.exon_sites
#PBS -o 16.exon_sites.o
#PBS -e 16.exon_sites.e
#PBS -P 11003054

# Define the directories and files
VCF_DIR="/fullgenotype/new_version/15.after_imputation_qc"
BED_FILE="/genotype/exon_sites/filtered/hg38exonsUCSC_modifiedCHR.bed"
OUTPUT_DIR="/project/adipose/fullgenotype/new_version/16.exon_sites"

#Loop through all chromosomes 1 to 22

for i in `seq 1 22`; do
    # Use the original VCF files directly
    ~/miniconda3/envs/qtl/bin/vcftools --gzvcf ${VCF_DIR}/chr${i}.dose_biallelic_filtered.vcf.gz --bed $BED_FILE --recode --recode-INFO-all --out ${OUTPUT_DIR}/chr${i}.dose_biallelic.filter.exon

    # Count the number of SNP variants and append to a file
    SNP_COUNT=$(~/miniconda3/envs/qtl/bin/bcftools view -v snps ${OUTPUT_DIR}/chr${i}.dose_biallelic.filter.exon.recode.vcf | grep -v '^#' | wc -l)
    echo -e "chr${i}\t${SNP_COUNT}" >> ${OUTPUT_DIR}/00.genotype_exon_num_per_chr.txt

done
