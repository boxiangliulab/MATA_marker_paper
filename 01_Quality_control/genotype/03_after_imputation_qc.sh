#!/bin/bash
#PBS -q normal
#PBS -l select=1:ncpus=128:mem=50G
#PBS -l walltime=3:00:00
#PBS -N after_imputation_qc
#PBS -o after_imputation_qc.o
#PBS -e after_imputation_qc.e
#PBS -P 11003054

set -euo pipefail
module load parallel

BCFTOOLS="$HOME/miniconda3/envs/qtl/bin/bcftools"
input_dir="/project/adipose/fullgenotype/new_version/14.imputation_results"
output_dir="/project/adipose/fullgenotype/new_version/15.after_imputation_qc"

# One sample ID per line; IDs must exactly match VCF sample IDs.
CHINESE_SAMPLES="/project/adipose/fullgenotype/new_version/sample_by_ancestry/chinese.samples.txt"
MALAY_SAMPLES="/project/adipose/fullgenotype/new_version/sample_by_ancestry/malay.samples.txt"
INDIAN_SAMPLES="/project/adipose/fullgenotype/new_version/sample_by_ancestry/indian.samples.txt"

mkdir -p "$output_dir"

for f in "$CHINESE_SAMPLES" "$MALAY_SAMPLES" "$INDIAN_SAMPLES"; do
    [[ -s "$f" ]] || { echo "ERROR: sample list not found or empty: $f" >&2; exit 1; }
done

process_file() {
    file="$1"
    output_dir="$2"
    BCFTOOLS="$3"
    CHINESE_SAMPLES="$4"
    MALAY_SAMPLES="$5"
    INDIAN_SAMPLES="$6"

    filename=$(basename "$file" .vcf.gz)
    workdir="$output_dir/${filename}.tmp"
    mkdir -p "$workdir"

    echo "[INFO] Processing $filename"

    # 1. Pooled/basic post-imputation QC only:
    #    biallelic SNPs, R2 >= 0.8, missingness <= 5%.
    "$BCFTOOLS" view -m2 -M2 -v snps "$file" -Ou | \
      "$BCFTOOLS" +fill-tags -Ou -- -t F_MISSING | \
      "$BCFTOOLS" view -i 'INFO/R2>=0.8 && INFO/F_MISSING<=0.05' \
        -Oz -o "$workdir/basic_qc.vcf.gz"
    "$BCFTOOLS" index -t "$workdir/basic_qc.vcf.gz"

    # 2. Per-ancestry MAF/HWE QC.
    for anc in Chinese Malay Indian; do
        case "$anc" in
            Chinese) samples="$CHINESE_SAMPLES" ;;
            Malay)   samples="$MALAY_SAMPLES" ;;
            Indian)  samples="$INDIAN_SAMPLES" ;;
        esac

        "$BCFTOOLS" view -S "$samples" "$workdir/basic_qc.vcf.gz" -Ou | \
          "$BCFTOOLS" +fill-tags -Ou -- -t MAF,HWE | \
          "$BCFTOOLS" view -i 'INFO/MAF>=0.05 && INFO/HWE>=1e-6' \
            -Oz -o "$workdir/${anc}.pass.vcf.gz"
        "$BCFTOOLS" index -t "$workdir/${anc}.pass.vcf.gz"

        "$BCFTOOLS" query -f '%CHROM\t%POS\t%REF\t%ALT\n' \
          "$workdir/${anc}.pass.vcf.gz" > "$workdir/${anc}.pass.sites.tsv"
    done

    # 3. UNION across ancestries:
    cat "$workdir/Chinese.pass.sites.tsv" \
        "$workdir/Malay.pass.sites.tsv" \
        "$workdir/Indian.pass.sites.tsv" | \
        sort -k1,1V -k2,2n -k3,3 -k4,4 -u > "$workdir/union.pass.sites.tsv"

    # 4. Apply the union site list back to the pooled Asian VCF,
    "$BCFTOOLS" view -T "$workdir/union.pass.sites.tsv" \
      "$workdir/basic_qc.vcf.gz" \
      -Oz -o "$output_dir/${filename}_biallelic_filtered.vcf.gz"
    "$BCFTOOLS" index -t "$output_dir/${filename}_biallelic_filtered.vcf.gz"

    n_basic=$("$BCFTOOLS" view -H "$workdir/basic_qc.vcf.gz" | wc -l)
    n_chinese=$(wc -l < "$workdir/Chinese.pass.sites.tsv")
    n_malay=$(wc -l < "$workdir/Malay.pass.sites.tsv")
    n_indian=$(wc -l < "$workdir/Indian.pass.sites.tsv")
    n_union=$(wc -l < "$workdir/union.pass.sites.tsv")
    n_final=$("$BCFTOOLS" view -H "$output_dir/${filename}_biallelic_filtered.vcf.gz" | wc -l)

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "$filename" "$n_basic" "$n_chinese" "$n_malay" "$n_indian" "$n_union" "$n_final" \
      > "$output_dir/${filename}.counts.tsv"

    rm -rf "$workdir"
    echo "[INFO] Finished $filename: $n_final variants retained"
}

export -f process_file

find "$input_dir" -name "*.dose.vcf.gz" -print0 | \
  parallel -0 process_file {} "$output_dir" "$BCFTOOLS" \
    "$CHINESE_SAMPLES" "$MALAY_SAMPLES" "$INDIAN_SAMPLES"

# 5. Count the final VCFs
summary="$output_dir/snp_count_summary.txt"
printf 'File\tBasic_QC\tChinese_PASS\tMalay_PASS\tIndian_PASS\tUnion_PASS\tFinal\n' > "$summary"
cat "$output_dir"/*.counts.tsv >> "$summary"

total=$(awk 'NR>1 {sum+=$7} END {print sum+0}' "$summary")
echo -e "TOTAL_FINAL_SNPS\t$total" >> "$summary"

echo "[INFO] QC complete. Final SNP total: $total"
echo "[INFO] Summary: $summary"
