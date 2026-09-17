#!/bin/bash
# Barcode lists for every meta-cell and every donor.
#
# Convert the meta-cell map of each cell type and depot into the per-sample
# barcode lists that the junction-counting step (step 06/07) consumes:
#
#   ${METACELL_DIR}/metadata/<cell_type>/<depot>/<library>_<donor>/<meta_cell>.txt
#   ${METACELL_DIR}/metadata/<cell_type>/all/<depot>/<library>_<donor>/<donor>.txt
#
# The first grouping gives the meta-cell level phenotypes used for sQTL
# mapping; the second pools all nuclei of a donor into one pseudobulk sample
# ("all"), whose summary statistics are used to filter the splice sites.

set -euo pipefail
cd "$(dirname "$0")"

# Paths
METACELL_DIR=../data/03_meta_cell

for depot in SAT VAT; do
    echo "=== ${depot} ==="

    cut -f1 "../reference/celltypes.QTL.${depot}.txt" | while read -r ct; do
        echo "--- ${ct} ---"

        out_dir="${METACELL_DIR}/metadata/${ct}"
        mkdir -p "${out_dir}"

        # The meta-cell size is part of the file name and differs per cell type.
        meta_cell_csv=$(ls "${METACELL_DIR}/KNN/${ct}/${ct}_${depot}_meta_cell_knn_fixed_"*"_centroid.csv")

        # One barcode list per meta-cell.
        python3 write_barcode_lists.py \
            --mode meta_cell \
            --input ${meta_cell_csv} \
            --out_dir "${out_dir}"

        # One barcode list per pseudobulk sample (all nuclei of that donor).
        python3 write_barcode_lists.py \
            --mode pseudobulk \
            --input ${meta_cell_csv} \
            --out_dir "${out_dir}"
    done
done

echo "Barcode lists written to ${METACELL_DIR}/metadata"
