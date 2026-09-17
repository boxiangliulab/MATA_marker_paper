#!/bin/bash
# Junction counts for every meta-cell.
#
# Count the junction reads of stage 02 over the nuclei of each meta-cell
# (ISSAC juncstat). Junction-spanning reads are grouped by genomic
# coordinates, strand, cell barcode and UMI, and each junction is supported by
# the number of unique CB-UMI combinations, so PCR duplicates are not counted
# twice. These counts are the input of the phenotype (PSI) construction in
# stage 04.

set -euo pipefail
cd "$(dirname "$0")"

# Paths
JUNC_DIR=../data/02_junction
METACELL_DIR=../data/03_meta_cell

for depot in SAT VAT; do
    echo "=== ${depot} ==="

    cut -f1 "../reference/celltypes.QTL.${depot}.txt" | while read -r ct; do
        echo "--- ${ct} ---"

        meta_dir="${METACELL_DIR}/metadata/${ct}/${depot}"
        out_dir="${METACELL_DIR}/juncstat/${ct}/${depot}"
        mkdir -p "${out_dir}"

        # One sub-directory per sample (<library>_<donor>) ...
        for sample_path in "${meta_dir}"/*/; do
            sample=$(basename "${sample_path}")

            # ... holding one barcode list per meta-cell.
            for barcode_file in "${sample_path}"*; do
                meta_cell=$(basename "${barcode_file}")
                ISSAC juncstat \
                    "${barcode_file}" \
                    "${JUNC_DIR}/${sample}.junc_barcoded" \
                    "${out_dir}/${sample}_${meta_cell}.stat"
            done
        done
    done
done

echo "Meta-cell junction counts written to ${METACELL_DIR}/juncstat"
