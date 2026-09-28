#!/bin/bash
#PBS -q normal
#PBS -l select=1:ncpus=2:mem=300G
#PBS -l walltime=2:00:00
#PBS -N 00.vcf_bygroup
#PBS -o 00.vcf_bygroup.o
#PBS -e 00.vcf_bygroup.e
#PBS -P 11003054

WORKDIR="/project/adipose/genotype/02.vcf_bygroup"
GROUP_DIR="/project/adipose/genotype/01.group_list"
VCF_FILE="../00.merge_genotype/adipose_after_imputation_300_merged_exon_genome_rename.vcf"

cd $WORKDIR

for group_file in ${GROUP_DIR}/*.txt
do
    group_name=$(basename "$group_file" .txt)
    output_file="${group_name}_GT.vcf"
    ~/miniconda3/envs/qtl/bin/bcftools view -S "${group_file}" -o - -O v "${VCF_FILE}" | \
    ~/miniconda3/envs/qtl/bin/bcftools annotate -x ^FORMAT/GT -o "${output_file}" -O v

    echo "processing: ${output_file}"
done

echo "done"
