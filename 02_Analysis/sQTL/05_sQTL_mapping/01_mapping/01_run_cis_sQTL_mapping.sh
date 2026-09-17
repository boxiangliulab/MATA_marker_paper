#!/bin/bash
# Cis-sQTL mapping.
#
# Map cis-sQTL for every splice site, cell type and depot with ISSAC.
#
# For each splice site ISSAC model fits a mixed-effects null model with
# fixed effects for the covariates and principal components and a random
# genetic effect given by the GRM, which accounts for relatedness and
# residual population structure. ISSAC QTL then score-tests every variant
# within +/-1 Mb of the splice site against that null model.
#
# Covariates (prepared in ../00_covariates): 5 splicing PCs, 10 genotype PCs
# in SAT and 6 in VAT, plus age, sex, diabetes status and BMI.
#
# Splice sites are processed chromosome by chromosome and, within a
# chromosome, in chunks of ${CHUNK_SIZE} sites; each chunk is handed to the
# worker script run_cis_sQTL_chunk.sh.

set -euo pipefail
cd "$(dirname "$0")"

# Paths
PHENO_DIR=../../data/04_phenotype
MODEL_DIR=../../data/05_sQTL_mapping

CHUNK_SIZE=100

log() { echo "[$(date '+%F %T')] $*"; }

for depot in SAT VAT; do
    log "DEPOT=${depot}"
    if [ "${depot}" == "SAT" ]; then g=10; else g=6; fi

    while IFS= read -r ct; do
        [ -z "${ct}" ] && continue
        log "  CELL TYPE=${ct}"

        pheno_file="${PHENO_DIR}/phenotype_output/${ct}/${depot}/phenotype_output.noFilter_rename.filtered"
        pc_file="${MODEL_DIR}/PCs/${ct}/${depot}/all_PCs.g${g}.p5.format.txt"
        model_dir="${MODEL_DIR}/model/${ct}/${depot}/p5_g${g}"
        out_dir="${MODEL_DIR}/cisQTL/${ct}/${depot}/p5_g${g}"
        chunk_dir="${MODEL_DIR}/model/${ct}/${depot}/chunks"

        [ -s "${pheno_file}" ] || { log "  [WARN] missing ${pheno_file}"; continue; }
        [ -s "${pc_file}" ]    || { log "  [WARN] missing ${pc_file}"; continue; }

        mkdir -p "${model_dir}" "${out_dir}" "${chunk_dir}"

        # Sample names shared by the phenotype and covariate files; the
        # covariate file header is one sample per column.
        head -1 "${pc_file}" | tr '\t' '\n' > "${out_dir}/sample_names.common"

        for chr in $(seq 1 22); do
            log "    CHR=${chr}"

            # Split the splice sites of this chromosome (the phenotype IDs start
            # with "chr<N>:") into chunks of ${CHUNK_SIZE} lines.
            chunk_prefix="${chunk_dir}/chunk_chr${chr}_"
            rm -f "${chunk_prefix}"* 2>/dev/null || true

            num_phenos=$(
                awk -v chr="chr${chr}:" -v cs="${CHUNK_SIZE}" -v pfx="${chunk_prefix}" '
                    NR == 1 { next }                       # skip the header
                    $0 ~ "^" chr {
                        if (NF == 0) next
                        c++
                        f = sprintf("%s%03d", pfx, int((c - 1) / cs))
                        print > f
                    }
                    END { print c + 0 }
                ' "${pheno_file}"
            )

            [ "${num_phenos}" -gt 0 ] || continue
            log "      splice sites: ${num_phenos}"

            for chunk_file in "${chunk_prefix}"*; do
                [ -s "${chunk_file}" ] || continue
                bash run_cis_sQTL_chunk.sh \
                    "${ct}" "${depot}" "${chr}" "${chunk_file}"
            done
        done
    done < <(cut -f1 "../../reference/celltypes.QTL.${depot}.txt")
done

log "cis-sQTL mapping finished: ${MODEL_DIR}/cisQTL"
