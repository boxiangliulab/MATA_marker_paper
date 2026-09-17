#!/bin/bash
# Genotype PC table.
#
# Extract the genotype principal components of each donor from the cohort
# PCA and write one table per number of PCs used. They account for population
# structure in the cis-sQTL model; the number was chosen per depot, 10 for
# SAT and 6 for VAT.

set -euo pipefail
cd "$(dirname "$0")"

# Paths
# Ancestry-specific genotype PCA of the cohort (col 3 donor ID, PCs from col 4)
GENOTYPE_PCA_FILE=../../data/genotype/genotype_PCA.csv
MODEL_DIR=../../data/05_sQTL_mapping

pc_dir="${MODEL_DIR}/PCs"
mkdir -p "${pc_dir}"

for n in 10 6; do
    out="${pc_dir}/genotype_PCs.${n}.txt"

    # Donor ID (column 3) plus the first n principal components (columns 4..n+3)
    last_col=$((n + 3))

    # Header: donor_id + gPC1..gPCn
    header="donor_id"
    for j in $(seq 1 "${n}"); do
        header="${header}\tgPC${j}"
    done
    echo -e "${header}" > "${out}"

    cut -d',' -f3-"${last_col}" "${GENOTYPE_PCA_FILE}" \
        | tail -n +2 \
        | tr ',' '\t' \
        | tr -d '"' >> "${out}"

    echo "Wrote ${out}"
done
