#!/bin/bash
#PBS -q normal
#PBS -l select=1:ncpus=24:mem=300G
#PBS -l walltime=2:00:00
#PBS -N beforeimpu
#PBS -o beforeimpu.o
#PBS -e beforeimpu.e
#PBS -P 11003054

# Before imputation Genotype QC pipeline (Suitable for GWAS, QTL analysis and demultiplex)
# Author: Wang Wenjing
# Version: 2.0

# =============== Initial Preparation ===============
# Path and Parameter should be defined firstly
output_dir="/project/adipose/fullgenotype/new_version"  # Working Path
ref_fasta="/share/data/reference/hg38/fa_ENSEMBL/Homo_sapiens.GRCh38.dna.primary_assembly.fa"  # hg38 reference genome
snpflip_exec=".conda/envs/python2.7/bin/snpflip"  # snpflip
plink_exec=".conda/envs/muscle/bin/plink" # plink
bcftools_exec=".conda/envs/muscle/bin/bcftools" # bcftools
liftover=".conda/envs/muscle/bin/liftOver" # LiftOver
chain_file="/data/projects/11003054/e1101919/software/hg19ToHg38.over.chain.gz" # hg19ToHg38 chain
tabix=".conda/envs/muscle/bin/tabix" # tabix
bgzip=".conda/envs/muscle/bin/bgzip" # bgzip
Rscript=".conda/envs/muscle/bin/Rscript" # Rscript
# =============== Step 2: Keep valid alleles ===============
echo "[INFO] 02.Keep valid alleles..."
valid_folder="$output_dir/02.valid_alleles"
mkdir -p $valid_folder

$plink_exec --bfile $output_dir/01.merge_300/merged_bfiles \
            --recode vcf-iid \
            --out $valid_folder/merged_bfiles

echo "[INFO] Conversion to VCF completed."

$bcftools_exec view -m2 -M2 -v snps $valid_folder/merged_bfiles.vcf \
             -Oz -o $valid_folder/filtered_snps.vcf.gz
$bcftools_exec index $valid_folder/filtered_snps.vcf.gz

$plink_exec --vcf $valid_folder/filtered_snps.vcf.gz \
            --make-bed \
            --const-fid \
            --out $valid_folder/valid_snps

valid_snp_count=$(wc -l < "$valid_folder/valid_snps.bim")

echo "[INFO] Step2. Keep valid alleles completed, total valid SNPs: $valid_snp_count"

# =============== Step 3: Sex check ===============
echo "[INFO] 03. Sex check..."
sex_folder="$output_dir/03.sex_check"
mkdir -p $sex_folder

if [[ -n "${reported_sex_file:-}" && -s "$reported_sex_file" ]]; then
    echo "[INFO] Updating sex information..."
    $plink_exec --bfile $valid_folder/valid_snps \
                --update-sex $reported_sex_file \
                --make-bed \
                --out $sex_folder/updated_sex

    input_file="$sex_folder/updated_sex"
else
    echo "[WARNING] No reported sex file provided. Skipping sex update..."
    input_file="$valid_folder/valid_snps"
fi

$plink_exec --bfile $input_file \
            --check-sex \
            --out $sex_folder/sex_check

awk 'NR>1 {print $1, $2, $3, $4, ($3==$4 ? "match" : "mismatch")}' $sex_folder/sex_check.sexcheck > $sex_folder/sex_comparison.txt

$Rscript /data/projects/11003054/e1101919/packed_GT_processed/plot_sex_check.R $sex_folder/sex_check.sexcheck $sex_folder/sex_check_plot.pdf

echo "[INFO] Step3. Sex check completed."

# =============== Step 4: Extract autosomes and bfiles2VCF ===============
echo "[INFO] 04.Extract autosomes and convert bfiles to VCF files..."
autosome_folder="$output_dir/04.Autosomes_VCF"
mkdir -p $autosome_folder
$plink_exec --bfile $valid_folder/valid_snps \
            --chr 1-22 \
            --make-bed \
            --out $autosome_folder/autosomes

auto_snp_count=$(wc -l < "$autosome_folder/autosomes.bim")

echo "[INFO] Step4. Extract autosomes and bfiles2VCF completed, total SNPs on autosome: $auto_snp_count"

# If ERRPR, check the "chr" perfix please.

# =============== Step 5: Convert VCF from hg19 to hg38 ===============
echo "[INFO] 05.Converting VCF files to hg38 format..."
hg38_folder="$output_dir/05.hg38_VCF"
mkdir -p $hg38_folder

$plink_exec --bfile $autosome_folder/autosomes \
            --recode tab \
            --const-fid \
            --out $hg38_folder/for_liftover

awk '{print "chr" $1, $4, $4+1, $2}' $hg38_folder/for_liftover.map > $hg38_folder/for_liftover.bed

$liftover $hg38_folder/for_liftover.bed $chain_file $hg38_folder/lifted_over.bed $hg38_folder/unmapped.bed
awk '{gsub("chr", "", $1); print $4, $1, int($3)}' $hg38_folder/lifted_over.bed > $hg38_folder/lifted_over.pos

# Remove unmapped SNPs
cut -f4 $hg38_folder/unmapped.bed > $hg38_folder/unmapped_ids.txt
$plink_exec --bfile $autosome_folder/autosomes \
            --exclude $hg38_folder/unmapped_ids.txt \
            --make-bed \
            --out $hg38_folder/filtered_snps

# update bp information
awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$4]=$2; next} {if ($2 in a) $4=a[$2]}1' $hg38_folder/lifted_over.bed $hg38_folder/filtered_snps.bim > $hg38_folder/autosomes_hg38_bfiles.bim

scp $hg38_folder/filtered_snps.bed $hg38_folder/autosomes_hg38_bfiles.bed
scp $hg38_folder/filtered_snps.fam $hg38_folder/autosomes_hg38_bfiles.fam

echo "[INFO] Step5. Convert VCF from hg19 to hg38 completed."

# =============== Step 6: Sample-level QC filtering ===============
echo "[INFO] Sample-level Quality Control..."
sample_qc_folder="$output_dir/06.sample_qc"
mkdir -p $sample_qc_folder

# Missingness check

$plink_exec --bfile $hg38_folder/autosomes_hg38_bfiles \
            --missing \
            --out $sample_qc_folder/missing

# Missingness filtering
$plink_exec --bfile $hg38_folder/autosomes_hg38_bfiles --mind 0.02 --make-bed --out $sample_qc_folder/missingness_filtered_data

# Missingness Plots
$Rscript /data/projects/11003054/e1101919/packed_GT_processed/plot_missing.R $sample_qc_folder/missing.imiss $sample_qc_folder/sample_missingness.pdf
$Rscript /data/projects/11003054/e1101919/packed_GT_processed/cumulative_snpmiss.R $sample_qc_folder/missing.imiss $sample_qc_folder/cumulative_sample_missingness.pdf

# Heterozygous ratio check
$plink_exec --bfile $sample_qc_folder/missingness_filtered_data \
            --het \
            --out $sample_qc_folder/het

# Het Plots
$Rscript /data/projects/11003054/e1101919/packed_GT_processed/heterozygosity_outliers_list.R $sample_qc_folder/het.het $sample_qc_folder/fail-het-qc.txt

# missingness VS. het plot
$Rscript /data/projects/11003054/e1101919/packed_GT_processed/missVShet.R $sample_qc_folder/missing.imiss $sample_qc_folder/het.het $sample_qc_folder/missingnessVShet.pdf
QCed_sample_count=$(wc -l < "$sample_qc_folder/missingness_filtered_data.fam")
echo "[INFO] Step6. Sample-level QC filtering completed, $QCed_sample_count samples pass the QC."

# =============== Step 7: SNP-level QC filtering ===============
echo "[INFO] SNP-level Quality Control..."
maf_hwe_folder="$output_dir/07.snp_qc"
mkdir -p $maf_hwe_folder

# Missingness filter
$Rscript /data/projects/11003054/e1101919/packed_GT_processed/plot_missing.R $sample_qc_folder/missing.lmiss $maf_hwe_folder/variant_missing.pdf
$Rscript /data/projects/11003054/e1101919/packed_GT_processed/cumulative_snpmiss.R $sample_qc_folder/missing.lmiss $maf_hwe_folder/cumulative_variant_missingness.pdf

$plink_exec --bfile $sample_qc_folder/missingness_filtered_data --geno 0.05 --make-bed --out $maf_hwe_folder/missingness_filtered

# MAF filter
$plink_exec --bfile $maf_hwe_folder/missingness_filtered --freq --out $maf_hwe_folder/maf_filtered

$Rscript /data/projects/11003054/e1101919/packed_GT_processed/plot_maf_density.R $maf_hwe_folder/maf_filtered.frq $maf_hwe_folder/maf_density.pdf

$plink_exec --bfile $maf_hwe_folder/missingness_filtered \
            --maf 0.01 \
            --make-bed \
            --out $maf_hwe_folder/maf_filtered

# HWE filter
$plink_exec --bfile $maf_hwe_folder/maf_filtered --hardy --out $maf_hwe_folder/hwe_filtered

$Rscript /data/projects/11003054/e1101919/packed_GT_processed/hwe_manhattan.R $maf_hwe_folder/hwe_filtered.hwe $maf_hwe_folder/maf_filtered.bim $maf_hwe_folder/hwe_manhattan.pdf

$plink_exec --bfile $maf_hwe_folder/maf_filtered \
            --hwe 1e-6 \
            --make-bed \
            --out $maf_hwe_folder/hwe_filtered

QCed_snp_count=$(wc -l < "$maf_hwe_folder/hwe_filtered.bim")

echo "[INFO] Step7. SNP-level QC filtering completed, $QCed_snp_count SNPs pass the QC."

# =============== Step 8: Relatedness Check ===============
echo "[INFO] 08. Checking relatedness among samples..."
relatedness_folder="$output_dir/08.relatedness"
mkdir -p $relatedness_folder

# Calculate the relatedness between samples
$plink_exec --bfile $maf_hwe_folder/hwe_filtered \
            --genome \
            --out $relatedness_folder/relatedness

# Samples with higher relatedness (pi_hat > 0.2)
awk '$10 > 0.2' $relatedness_folder/relatedness.genome > $relatedness_folder/close_relatives.txt

# Relatedness plot
$Rscript /data/projects/11003054/e1101919/packed_GT_processed/plot_relatedness.R $relatedness_folder/relatedness.genome $relatedness_folder/close_relatives.txt $relatedness_folder/relatedness_plot.pdf

echo "[INFO] Step8. Relatedness Check completed."

# =============== Step 9: SNPflip ===============
echo "[INFO] 09. Running SNPflip..."
snpflip_folder="$output_dir/09.snpflip"
mkdir -p $snpflip_folder

$snpflip_exec --bim $maf_hwe_folder/hwe_filtered.bim \
              --fasta $ref_fasta \
              --output $snpflip_folder/snps_flip

# Flip the reverse SNPs and delete the ambiguous SNPs
$plink_exec --bfile $maf_hwe_folder/hwe_filtered \
            --flip $snpflip_folder/snps_flip.reverse \
            --exclude $snpflip_folder/snps_flip.ambiguous \
            --make-bed \
            --out $snpflip_folder/filtered_forward_snps

echo "[INFO] Step9. SNPflip completed."

# =============== Step 10: Convert bfiles to VCF ===============
echo "[INFO] 10. Converting bfiles to VCF..."
forward_vcf_folder="$output_dir/10.forward_vcf"
mkdir -p $forward_vcf_folder

$plink_exec --bfile $snpflip_folder/filtered_forward_snps \
            --recode vcf \
            --out $forward_vcf_folder/forward_snps

# BGzip compress the VCF files
$bgzip -c $forward_vcf_folder/forward_snps.vcf > $forward_vcf_folder/forward_snps.vcf.gz
$tabix -p vcf $forward_vcf_folder/forward_snps.vcf.gz
echo "[INFO] Step10. Convert bfiles to VCF completed."

# =============== Step 11: Normalize Alleles using bcftools ===============
echo "[INFO] 11. Normalizing alleles in VCF using bcftools..."
normalized_vcf_folder="$output_dir/11.normalized_vcf"
mkdir -p $normalized_vcf_folder

$bcftools_exec norm -f $ref_fasta \
                    -c ws \
                    -m -any \
                    -Oz -o $normalized_vcf_folder/normalized.vcf.gz \
                    $forward_vcf_folder/forward_snps.vcf.gz
echo "[INFO] Step11. Normalize Alleles using bcftools completed."

# =============== Step 12: Before Imputation QC ===============
echo "[INFO] 12. Performing QC before imputation..."
imputation_qc_folder="$output_dir/12.imputation_qc"
mkdir -p $imputation_qc_folder

$plink_exec --bfile $snpflip_folder/filtered_forward_snps \
            --missing \
            --out $imputation_qc_folder/missingness

$plink_exec --bfile $snpflip_folder/filtered_forward_snps \
            --freq \
            --out $imputation_qc_folder/maf

$plink_exec --bfile $snpflip_folder/filtered_forward_snps \
            --hardy \
            --out $imputation_qc_folder/hwe

echo "[INFO] Step12. Before Imputation QC completed."

# =============== Step 13: Prepare VCF for Imputation ===============
echo "[INFO] 13. Preparing VCF for imputation..."
imputation_vcf_folder="$output_dir/13.imputation_vcf"
mkdir -p $imputation_vcf_folder

# Ensure the "chr" prefix for the VCF files
zcat $normalized_vcf_folder/normalized.vcf.gz | sed '/^#/!s/^\([0-9XYM]\)/chr\1/' > $imputation_vcf_folder/chr_normalized.vcf
$bgzip -c $imputation_vcf_folder/chr_normalized.vcf > $imputation_vcf_folder/chr_normalized.vcf.gz
$tabix -p vcf $imputation_vcf_folder/chr_normalized.vcf.gz

# Split files into chr1~chr22
for chrom in {1..22}; do
    echo "[INFO] Splitting chromosome $chrom..."
    $bcftools_exec view -r chr$chrom $imputation_vcf_folder/chr_normalized.vcf.gz \
                        -Oz -o $imputation_vcf_folder/chr$chrom.vcf.gz
    $tabix -p vcf $imputation_vcf_folder/chr$chrom.vcf.gz
done
echo "[INFO] 13. Preparing VCF for imputation completed"
# =============== Completion ===============
# The files in $imputation_vcf_folder are ready for Michigan Imputation Server!
echo "[INFO] Congratulations! All finished! Results saved in  $output_dir"
