#!/bin/bash
# Aggregate results, call sQTL.
#
# Concatenate the per-splice-site association files produced in step 01 into
# one file per cell type and depot, then apply the two-level multiple-testing
# correction implemented in multiple_testing_correction.R:
#   Bonferroni within each splice site -> lead variant,
#   Benjamini-Hochberg across splice sites -> significant sQTL (FDR < 0.05).

set -euo pipefail
cd "$(dirname "$0")"

# Paths
MODEL_DIR=../../data/05_sQTL_mapping

results_dir="${MODEL_DIR}/results"
mkdir -p "${results_dir}"

for depot in SAT VAT; do
    echo "=== ${depot} ==="
    if [ "${depot}" == "SAT" ]; then g=10; else g=6; fi

    cut -f1 "../../reference/celltypes.QTL.${depot}.txt" | while read -r ct; do
        echo "--- ${ct} ---"

        in_dir="${MODEL_DIR}/cisQTL/${ct}/${depot}/p5_g${g}"
        merged_file="${results_dir}/${ct}.${depot}.p5_g${g}.result"

        n_files=$(find "${in_dir}" -maxdepth 1 -type f -name "*result" | wc -l)
        if [ "${n_files}" -eq 0 ]; then
            echo "  Warning: no *result files in ${in_dir}, skipping"
            continue
        fi

        # Concatenate in a reproducible order (find | sort handles the large
        # number of files that would overflow the shell's argument list).
        : > "${merged_file}"
        find "${in_dir}" -maxdepth 1 -type f -name "*result" -print0 \
            | sort -z \
            | xargs -0 cat >> "${merged_file}"

        echo "  Bonferroni correction within splice site + global FDR..."
        Rscript multiple_testing_correction.R \
            "${merged_file}" \
            "${merged_file}.lead_variants.fdr" \
            "${merged_file}.sQTL.fdr05"
    done
done

echo "Aggregated results written to ${results_dir}"
