#!/bin/bash
# Index the per-donor alignments.
#
# Build BAM indexes for every donor-level STARsolo alignment; ISSAC junctools
# (step 02) requires indexed, coordinate-sorted BAM files.

set -euo pipefail
cd "$(dirname "$0")"

# Paths
META_DIR=../data/00_metadata
STARSOLO_DIR=../data/01_STARsolo

while read -r sample; do
    echo "=== ${sample} ==="
    samtools index -@ 24 \
        "${STARSOLO_DIR}/${sample}/${sample}_Aligned.sortedByCoord.out.bam"
done < "${META_DIR}/list_samples.txt"

echo "All BAM files indexed."
