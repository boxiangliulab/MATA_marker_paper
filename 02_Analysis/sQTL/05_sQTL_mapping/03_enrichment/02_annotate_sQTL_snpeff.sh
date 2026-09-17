#!/bin/bash
# SnpEff annotation of the sQTL.
#
# Annotate the significant sQTL variants with SnpEff, which adds the
# functional consequence of each variant (ANN= field) used by the enrichment
# analysis.

set -euo pipefail
cd "$(dirname "$0")"

# Paths
MODEL_DIR=../../data/05_sQTL_mapping
SNPEFF_JAR=../../data/tools/snpEff/snpEff.jar

snpeff_dir="${MODEL_DIR}/SNPeff"

if [ ! -f "${snpeff_dir}/sQTL_variants_FDR05.vcf" ]; then
    echo "ERROR: input VCF not found. Run 01_convert_sQTL_to_VCF.R first." >&2
    exit 1
fi

java -Xmx8g -jar "${SNPEFF_JAR}" \
    -v \
    -stats "${snpeff_dir}/snpEff_summary.html" \
    GRCh38.86 \
    "${snpeff_dir}/sQTL_variants_FDR05.vcf" \
    > "${snpeff_dir}/sQTL_variants_FDR05_annotated.vcf"

echo "Annotated VCF written to ${snpeff_dir}/sQTL_variants_FDR05_annotated.vcf"
