#!/bin/bash
# Group junctions into splice sites (ISSAC pheno_group).
#
# Define the splice-site partner relationships with ISSAC pheno_group: for
# each splice site, all junctions sharing the same donor or acceptor
# coordinate are grouped into competing splice site partners, which is what
# allows the percent-spliced-in (PSI) value of each splice site to be
# computed in step 03.
#
# Junctions supported by fewer than READ_THRESHOLD reads across the data set
# are dropped.
#
# The step is run twice per cell type and depot:
#   - meta-cell level : the samples used for sQTL mapping
#   - pseudobulk      : one sample per donor, used in step 03 to filter the
#                       splice sites

set -euo pipefail
cd "$(dirname "$0")"

# Paths
METACELL_DIR=../data/03_meta_cell
PHENO_DIR=../data/04_phenotype

READ_THRESHOLD=30

out_root="${PHENO_DIR}/phenotype_group"

for depot in SAT VAT; do
    echo "=== ${depot} ==="

    cut -f1 "../reference/celltypes.QTL.${depot}.txt" | while read -r ct; do
        echo "--- ${ct} ---"

        # --- meta-cell level ---------------------------------------------------
        out_dir="${out_root}/${ct}/${depot}"
        mkdir -p "${out_dir}"
        ISSAC pheno_group \
            "${PHENO_DIR}/input_sample_file/${ct}_${depot}.txt" \
            "${METACELL_DIR}/juncstat/${ct}/${depot}" \
            "${out_dir}/phenotype_group" \
            "${READ_THRESHOLD}" \
            "${out_dir}/out"

        # --- pseudobulk level -------------------------------------------------
        out_dir_all="${out_root}/${ct}/all/${depot}"
        mkdir -p "${out_dir_all}"
        ISSAC pheno_group \
            "${PHENO_DIR}/input_sample_file/${ct}_all_${depot}.txt" \
            "${METACELL_DIR}/juncstat/${ct}/all/${depot}" \
            "${out_dir_all}/phenotype_group" \
            "${READ_THRESHOLD}" \
            "${out_dir_all}/out"
    done
done

echo "Splice-site groups written to ${out_root}"
