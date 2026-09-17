#!/bin/bash
# Run the meta-cell assignment for every cell type.
#
# Driver for steps 01 and 02: for each depot and each cell type tested in
# that depot it
#   (a) summarises the number of nuclei per donor (step 01) to obtain the
#       target meta-cell size k = ceil(median nuclei per donor / 4), and
#   (b) runs the KNN meta-cell assignment (step 02) with that k, requiring
#       each meta-cell to contain at least 80% of k nuclei.

set -euo pipefail
cd "$(dirname "$0")"

# Paths
# scVI-integrated object per cell type: <cell_type>/<cell_type>.reintegrated.h5ad
SCRNA_OBJ_DIR=../data/snRNA/by_celltype
METACELL_DIR=../data/03_meta_cell

for depot in SAT VAT; do
    echo "=== ${depot} ==="

    cut -f1 "../reference/celltypes.QTL.${depot}.txt" | while read -r ct; do
        echo "--- ${ct} ---"

        h5ad="${SCRNA_OBJ_DIR}/${ct}/${ct}.reintegrated.h5ad"
        out_dir="${METACELL_DIR}/KNN/${ct}"
        median_table="${out_dir}/${ct}_median_cells_per_donor.tsv"
        mkdir -p "${out_dir}"

        # (a) Nuclei per donor and the derived meta-cell size.
        if [ ! -s "${median_table}" ]; then
            python3 01_summarize_cells_per_donor.py \
                --input_h5ad "${h5ad}" \
                --cell_type "${ct}" \
                --output_dir "${out_dir}"
        fi

        # Column 3 of the median table is ceil(median nuclei per donor / 4).
        k=$(grep -w "${depot}" "${median_table}" | cut -f3)
        # Minimum meta-cell size = 80% of k (integer ceiling).
        min_size=$(( (k * 80 + 99) / 100 ))

        echo "    target meta-cell size: ${k} (minimum ${min_size})"

        # (b) KNN meta-cell assignment.
        python3 02_assign_meta_cells_knn.py \
            --input_h5ad "${h5ad}" \
            --cell_type "${ct}" \
            --depot "${depot}" \
            --meta_size "${k}" \
            --min_final_size "${min_size}" \
            --output_dir "${out_dir}"
    done
done

echo "Meta-cell assignment complete: ${METACELL_DIR}/KNN"
