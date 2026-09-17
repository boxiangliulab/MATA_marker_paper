#!/bin/bash
# Map one chunk of splice sites (called by 01_run_cis_sQTL_mapping.sh).

set -euo pipefail
cd "$(dirname "$0")"

# Paths
PHENO_DIR=../../data/04_phenotype
MODEL_DIR=../../data/05_sQTL_mapping
GENOTYPE_BY_CHR_DIR=../../data/genotype/by_chr

ct="$1"
depot="$2"
chr="$3"
chunk_file="$4"

if [ "${depot}" == "SAT" ]; then g=10; else g=6; fi

pheno_file="${PHENO_DIR}/phenotype_output/${ct}/${depot}/phenotype_output.noFilter_rename.filtered"
pc_file="${MODEL_DIR}/PCs/${ct}/${depot}/all_PCs.g${g}.p5.format.txt"
model_dir="${MODEL_DIR}/model/${ct}/${depot}/p5_g${g}"
out_dir="${MODEL_DIR}/cisQTL/${ct}/${depot}/p5_g${g}"

# Scratch space for the single-splice-site input files.
work=$(mktemp -d "${TMPDIR:-/tmp}/issac_$(basename "${chunk_file}")_XXXXXX")
trap 'rm -rf "${work}"' EXIT

header="${work}/header.txt"
tmp_pheno="${work}/pheno.txt"
tmp_site="${work}/site.txt"

head -1 "${pheno_file}" > "${header}"

while IFS= read -r line || [ -n "${line}" ]; do
    [ -z "${line}" ] && continue

    pheno="${line%% *}"                          # splice-site ID

    # One-splice-site phenotype file and the matching site list.
    cat "${header}" > "${tmp_pheno}"
    printf '%s\n' "${line}" >> "${tmp_pheno}"
    printf '%s\n' "${pheno}" > "${tmp_site}"

    # --- Fit the null model for this splice site ----------------------------
    ISSAC model \
        -s "${tmp_pheno}" \
        -p "${pc_file}" \
        -g "${MODEL_DIR}/GRM.296.txt" \
        -n "296" \
        -u "${model_dir}" \
        -o "model.${pheno//[^a-zA-Z0-9_-]/_}"   # file-system safe name

    # --- Score-test all cis variants (+/-1 Mb) ------------------------------
    ISSAC QTL \
        -s "${tmp_site}" \
        -p "${model_dir}" \
        -m "${out_dir}/sample_names.common" \
        -c "chr${chr}" \
        -x "${pc_file}" \
        -v "${GENOTYPE_BY_CHR_DIR}/adipose.chr${chr}.bcf" \
        -w 1000000 \
        -o "${out_dir}" \
        -t 1
done < "${chunk_file}"
