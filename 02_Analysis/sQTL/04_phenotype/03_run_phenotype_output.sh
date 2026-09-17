#!/bin/bash
# Compute PSI phenotypes (ISSAC pheno_output).
#
# Turn the splice-site read counts of step 02 into percent-spliced-in (PSI)
# phenotype matrices with ISSAC pheno_output, which takes two filters:
#   SD      minimum standard deviation of PSI across samples
#   NA_PROP maximum fraction of samples with a missing value
#
# Two matrices are produced per cell type and depot:
#   - pseudobulk (SD = 0.01, NA_PROP = 0.4)
#       Splice site definition is made stringent on the pseudobulk summary
#       statistics: a site is retained if it varies across samples
#       (SD >= 0.01) and is non-missing in at least 60% of the samples.
#   - meta-cell level (SD = 0, NA_PROP = 1, i.e. no filtering)
#       All splice sites present in the aggregated data are retained here;
#       the matrix is restricted to the pseudobulk-filtered sites in step 04.

set -euo pipefail
cd "$(dirname "$0")"

# Paths
PHENO_DIR=../data/04_phenotype

group_root="${PHENO_DIR}/phenotype_group"
out_root="${PHENO_DIR}/phenotype_output"

for depot in SAT VAT; do
    echo "=== ${depot} ==="

    cut -f1 "../reference/celltypes.QTL.${depot}.txt" | while read -r ct; do
        echo "--- ${ct} ---"

        out_dir="${out_root}/${ct}/${depot}"
        mkdir -p "${out_dir}"

        # --- pseudobulk: filtered, defines the splice sites to test -----------
        ISSAC pheno_output \
            "${group_root}/${ct}/all/${depot}/phenotype_group.site" \
            "${out_dir}/phenotype_output.all.filtered" \
            "${out_dir}/phenotype_output.all.prop" \
            0.01 \
            0.4

        # --- meta-cell level: unfiltered --------------------------------------
        ISSAC pheno_output \
            "${group_root}/${ct}/${depot}/phenotype_group.site" \
            "${out_dir}/phenotype_output.noFilter" \
            "${out_dir}/phenotype_output.noFilter.prop" \
            0 \
            1
    done
done

echo "PSI phenotype matrices written to ${out_root}"
