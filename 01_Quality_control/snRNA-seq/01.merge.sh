#!/bin/bash
#PBS -q normal
#PBS -l select=1:ncpus=4:mem=300G
#PBS -l walltime=2:00:00
#PBS -N 00.merge
#PBS -o 00.merge.o
#PBS -e 00.merge.e
#PBS -P 11003054

WORKDIR="/project/adipose/genotype/00.merge_genotype"
cd $WORKDIR
~/miniconda3/envs/qtl/bin/bcftools concat -Oz -o adipose_after_imputation_300_merged_exon_genome.vcf.gz \
    /project/adipose/fullgenotype/new_version/16.exon_sites/chr1.dose_biallelic.filter.exon.recode.vcf \
    /project/adipose/fullgenotype/new_version/16.exon_sites/chr10.dose_biallelic.filter.exon.recode.vcf \
    /project/adipose/fullgenotype/new_version/16.exon_sites/chr11.dose_biallelic.filter.exon.recode.vcf \
    /project/adipose/fullgenotype/new_version/16.exon_sites/chr12.dose_biallelic.filter.exon.recode.vcf \
    /project/adipose/fullgenotype/new_version/16.exon_sites/chr13.dose_biallelic.filter.exon.recode.vcf \
    /project/adipose/fullgenotype/new_version/16.exon_sites/chr14.dose_biallelic.filter.exon.recode.vcf \
    /project/adipose/fullgenotype/new_version/16.exon_sites/chr15.dose_biallelic.filter.exon.recode.vcf \
    /project/adipose/fullgenotype/new_version/16.exon_sites/chr16.dose_biallelic.filter.exon.recode.vcf \
    /project/adipose/fullgenotype/new_version/16.exon_sites/chr17.dose_biallelic.filter.exon.recode.vcf \
    /project/adipose/fullgenotype/new_version/16.exon_sites/chr18.dose_biallelic.filter.exon.recode.vcf \
    /project/adipose/fullgenotype/new_version/16.exon_sites/chr19.dose_biallelic.filter.exon.recode.vcf \
    /project/adipose/fullgenotype/new_version/16.exon_sites/chr2.dose_biallelic.filter.exon.recode.vcf \
    /project/adipose/fullgenotype/new_version/16.exon_sites/chr20.dose_biallelic.filter.exon.recode.vcf \
    /project/adipose/fullgenotype/new_version/16.exon_sites/chr21.dose_biallelic.filter.exon.recode.vcf \
    /project/adipose/fullgenotype/new_version/16.exon_sites/chr22.dose_biallelic.filter.exon.recode.vcf \
    /project/adipose/fullgenotype/new_version/16.exon_sites/chr3.dose_biallelic.filter.exon.recode.vcf \
    /project/adipose/fullgenotype/new_version/16.exon_sites/chr4.dose_biallelic.filter.exon.recode.vcf \
    /project/adipose/fullgenotype/new_version/16.exon_sites/chr5.dose_biallelic.filter.exon.recode.vcf \
    /project/adipose/fullgenotype/new_version/16.exon_sites/chr6.dose_biallelic.filter.exon.recode.vcf \
    /project/adipose/fullgenotype/new_version/16.exon_sites/chr7.dose_biallelic.filter.exon.recode.vcf \
    /project/adipose/fullgenotype/new_version/16.exon_sites/chr8.dose_biallelic.filter.exon.recode.vcf \
    /project/adipose/fullgenotype/new_version/16.exon_sites/chr9.dose_biallelic.filter.exon.recode.vcf
 snp_count=$(~/miniconda3/envs/qtl/bin/bcftools view -v snps adipose_after_imputation_300_merged_exon_genome.vcf.gz | grep -v '^#' | wc -l)
 echo "$snp_count SNPs finaled." >> 00.genotype_num_merge.txt
 gunzip adipose_after_imputation_300_merged_exon_genome.vcf.gz
~/miniconda3/envs/qtl/bin/bcftools query -l adipose_after_imputation_300_merged_exon_genome.vcf  > samples.txt
