#!/bin/bash
# Build the library and donor manifests.
#
# The multiplexed snRNA-seq libraries were processed with Cell Ranger, one
# BAM per library. This script builds the two manifests that drive the rest
# of the demultiplexing stage:
#
#   1) ${META_DIR}/sample_list_group.tsv
#      <sequencing group> <TAB> <library ID>, one row per library.
#
#   2) ${META_DIR}/<library>.tsv
#      The subset of the cell-level metadata table belonging to that library
#      (with the original header kept), used downstream to pull each donor's
#      barcodes out of the pooled BAM.

set -euo pipefail
cd "$(dirname "$0")"

# Paths
# Cell Ranger output: <group>/<library>/outs/possorted_genome_bam.bam
CELLRANGER_DIR=../data/cellranger
# Cell-level metadata after library QC (col 1 barcode, 4 doublet score,
# 5 donor, 6 library)
CELL_METADATA=../data/snRNA/metadata_after_libraryQC.tsv
META_DIR=../data/00_metadata

mkdir -p "${META_DIR}"

sample_list="${META_DIR}/sample_list_group.tsv"

# --- 1. One row per library: sequencing group + library ID -------------------
# The group and library names are taken from the Cell Ranger directory layout
# (<group>/<library>/outs/possorted_genome_bam.bam).
# Pilot and add-on runs are excluded; only the full-scale libraries are used.
ls "${CELLRANGER_DIR}"/*/*/outs/possorted_genome_bam.bam | while read -r bam; do
    library=$(basename "$(dirname "$(dirname "${bam}")")")
    group=$(basename "$(dirname "$(dirname "$(dirname "${bam}")")")")
    printf '%s\t%s\n' "${group}" "${library}"
done | sort -k1,1V -k2,2V | grep -v pilot | grep -v addon > "${sample_list}"

echo "Libraries found: $(wc -l < "${sample_list}")"

# --- 2. Split the cell-level metadata table by library -----------------------
# Column 6 of ${CELL_METADATA} holds the library ID. The header line is kept so
# that the per-library files are self-describing.
cut -f2 "${sample_list}" | while read -r library; do
    {
        head -1 "${CELL_METADATA}"
        awk -v OFS='\t' -v lib="${library}" '$6 == lib' "${CELL_METADATA}"
    } > "${META_DIR}/${library}.tsv"
done

echo "Per-library metadata written to ${META_DIR}"
