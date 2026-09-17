#!/bin/bash
# Junction counts per pseudobulk sample.
#
# Same as step 05 (identical CB-UMI collapsed counting), but pooling all
# nuclei of a donor within a cell type and depot into a single pseudobulk
# sample instead of splitting them into meta-cells. The pseudobulk samples
# are deeper, and their summary statistics are used in stage 04 to filter the
# splice sites; the meta-cell counts are what the sQTL mapping uses.

set -euo pipefail
cd "$(dirname "$0")"

# Paths
JUNC_DIR=../data/02_junction
METACELL_DIR=../data/03_meta_cell

for depot in SAT VAT; do
    echo "=== ${depot} ==="

    cut -f1 "../reference/celltypes.QTL.${depot}.txt" | while read -r ct; do
        echo "--- ${ct} ---"

        meta_dir="${METACELL_DIR}/metadata/${ct}/all/${depot}"
        out_dir="${METACELL_DIR}/juncstat/${ct}/all/${depot}"
        mkdir -p "${out_dir}"

        for sample_path in "${meta_dir}"/*/; do
            sample=$(basename "${sample_path}")

            for barcode_file in "${sample_path}"*; do
                donor=$(basename "${barcode_file}")
                ISSAC juncstat \
                    "${barcode_file}" \
                    "${JUNC_DIR}/${sample}.junc_barcoded" \
                    "${out_dir}/${sample}_${donor}.stat"
            done
        done
    done
done

echo "Pseudobulk junction counts written to ${METACELL_DIR}/juncstat"
