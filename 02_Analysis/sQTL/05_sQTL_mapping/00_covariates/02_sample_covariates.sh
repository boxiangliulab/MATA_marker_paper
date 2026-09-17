#!/bin/bash
# Donor covariate table.
#
# Collect the donor-level covariates used in the cis-sQTL model -- age, sex,
# diabetes status and BMI -- from the frozen snRNA-seq metadata table, one
# row per donor.

set -euo pipefail
cd "$(dirname "$0")"

# Paths
# Frozen snRNA-seq metadata (cols 12 donor, 17 age, 18 sex, 20 diabetes, 21 BMI)
DONOR_METADATA=../../data/snRNA/metadata_freeze.tsv
MODEL_DIR=../../data/05_sQTL_mapping

pc_dir="${MODEL_DIR}/PCs"
mkdir -p "${pc_dir}"

out="${pc_dir}/sample_covariates.txt"

echo -e "donor_id\tAge\tSex\tDiabetic\tBMI" > "${out}"

# One row per nucleus in the input, so the donor rows are de-duplicated.
tail -n +2 "${DONOR_METADATA}" | cut -f 12,17,18,20,21 | sort -u >> "${out}"

echo "Wrote ${out} ($(($(wc -l < "${out}") - 1)) donors)"
