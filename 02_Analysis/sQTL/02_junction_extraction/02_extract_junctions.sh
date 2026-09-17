#!/bin/bash
# Barcode-aware splice-junction extraction.
#
# Extract split reads (spliced alignments) from each donor BAM with
# ISSAC junctools, keeping the cell barcode of every junction-spanning read.
# The resulting "barcoded junction" files are the raw material for the
# meta-cell junction counts computed in stage 03.
#
# ISSAC junctools parameters:
#   -a 8       minimum anchor length on either side of the junction
#   -m 50      minimum intron length
#   -M 500000  maximum intron length
#   -s FR      library strandedness

set -euo pipefail
cd "$(dirname "$0")"

# Paths
META_DIR=../data/00_metadata
STARSOLO_DIR=../data/01_STARsolo
JUNC_DIR=../data/02_junction

mkdir -p "${JUNC_DIR}"

while read -r sample; do
    echo "=== ${sample} ==="
    ISSAC junctools extract \
        -a 8 -m 50 -M 500000 \
        -s FR "${STARSOLO_DIR}/${sample}/${sample}_Aligned.sortedByCoord.out.bam" \
        -o "${JUNC_DIR}/${sample}.junc_barcoded"
done < "${META_DIR}/list_samples.txt"

echo "Barcoded junction files written to ${JUNC_DIR}"
