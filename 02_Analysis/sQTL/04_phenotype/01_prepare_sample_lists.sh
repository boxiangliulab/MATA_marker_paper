#!/bin/bash
# Sample lists for the phenotype construction.
#
# ISSAC pheno_group needs an explicit list of the samples to combine. For
# every cell type and depot this script lists the junction-count files
# produced in stage 03 and strips the .stat suffix, once for the meta-cell
# level samples and once for the pseudobulk samples.

set -euo pipefail
cd "$(dirname "$0")"

# Paths
METACELL_DIR=../data/03_meta_cell
PHENO_DIR=../data/04_phenotype

sample_file_dir="${PHENO_DIR}/input_sample_file"
mkdir -p "${sample_file_dir}"

for depot in SAT VAT; do
    echo "=== ${depot} ==="

    cut -f1 "../reference/celltypes.QTL.${depot}.txt" | while read -r ct; do
        echo "--- ${ct} ---"

        # Meta-cell level samples
        ls "${METACELL_DIR}/juncstat/${ct}/${depot}/" \
            | sed 's/\.stat$//' > "${sample_file_dir}/${ct}_${depot}.txt"

        # Pseudobulk samples
        ls "${METACELL_DIR}/juncstat/${ct}/all/${depot}/" \
            | sed 's/\.stat$//' > "${sample_file_dir}/${ct}_all_${depot}.txt"
    done
done

echo "Sample lists written to ${sample_file_dir}"
