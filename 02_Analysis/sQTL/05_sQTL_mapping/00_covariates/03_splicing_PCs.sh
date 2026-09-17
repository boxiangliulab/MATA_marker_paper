#!/bin/bash
# Splicing PCs per cell type.
#
# Run the PCA of the meta-cell PSI matrix for every cell type and depot and
# keep the top 5 components as covariates for the cis-sQTL model; they
# capture global variation in splicing patterns.

set -euo pipefail
cd "$(dirname "$0")"

# Paths
PHENO_DIR=../../data/04_phenotype
MODEL_DIR=../../data/05_sQTL_mapping

for depot in SAT VAT; do
    echo "=== ${depot} ==="

    cut -f1 "../../reference/celltypes.QTL.${depot}.txt" | while read -r ct; do
        pheno_file="${PHENO_DIR}/phenotype_output/${ct}/${depot}/phenotype_output.noFilter_rename.prop"

        if [ ! -s "${pheno_file}" ]; then
            echo "  skip ${ct}: ${pheno_file} not found"
            continue
        fi

        echo "--- ${ct} ---"
        out_dir="${MODEL_DIR}/PCs/${ct}/${depot}"
        mkdir -p "${out_dir}"

        python3 compute_splicing_PCs.py \
            --input "${pheno_file}" \
            --output "${out_dir}/splicing_PCs.5.txt" \
            --n_components "5"
    done
done

echo "Splicing PCs written under ${MODEL_DIR}/PCs"
