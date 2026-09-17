#!/bin/bash
# SnpEff annotation of the background.
#
# Annotate the background variants created in step 04 with SnpEff, using the
# same database as the significant sQTL variants, so that the two sets can be
# compared annotation type by annotation type in step 06.

set -euo pipefail
cd "$(dirname "$0")"

# Paths
MODEL_DIR=../../data/05_sQTL_mapping
SNPEFF_JAR=../../data/tools/snpEff/snpEff.jar

qtlbhm_dir="${MODEL_DIR}/qtlBHM"

if [ ! -f "${qtlbhm_dir}/background_variants.vcf" ]; then
    echo "ERROR: input VCF not found. Run 04_create_background_set.R first." >&2
    exit 1
fi

java -Xmx8g -jar "${SNPEFF_JAR}" \
    -v \
    -stats "${qtlbhm_dir}/background_snpEff_summary.html" \
    GRCh38.86 \
    "${qtlbhm_dir}/background_variants.vcf" \
    > "${qtlbhm_dir}/background_variants_annotated.vcf"

echo "Annotated VCF written to ${qtlbhm_dir}/background_variants_annotated.vcf"
