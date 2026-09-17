#!/bin/bash
# Split pooled libraries into per-donor FASTQ files.
#
# Each library is a pool of several donors. For every donor in every library
# we (i) extract that donor's nuclei from the pooled Cell Ranger BAM with
# subset-bam, and (ii) convert the resulting BAM back to FASTQ with
# bamtofastq, so that the reads can be re-aligned per donor (step 01 of
# 01_alignment) with donor-specific genotypes for WASP filtering.

set -euo pipefail
cd "$(dirname "$0")"

# Paths
CELLRANGER_DIR=../data/cellranger
META_DIR=../data/00_metadata
FASTQ_DIR=../data/00_bamtofastq

donor_list_dir="${META_DIR}/donor_list"
barcode_dir="${META_DIR}/barcodes"

mkdir -p "${donor_list_dir}" "${barcode_dir}" "${FASTQ_DIR}"

# Nuclei with a doublet score >= this value are discarded (column 4 of the
# cell-level metadata table).
DOUBLET_SCORE_MAX=0.85

while IFS=$'\t' read -r group library; do
    echo "=== ${group} / ${library} ==="

    meta_file="${META_DIR}/${library}.tsv"

    # --- Donors multiplexed in this library (column 5) -----------------------
    tail -n +2 "${meta_file}" \
        | awk -v max="${DOUBLET_SCORE_MAX}" '$4 < max' \
        | cut -f5 | sort -u > "${donor_list_dir}/${library}.txt"

    while read -r donor; do
        echo "--- donor ${donor} ---"
        full_name="${library}_${donor}"
        barcode_file="${barcode_dir}/${full_name}.txt"

        # --- Barcodes of this donor -------------------------------------------
        # Barcodes in the metadata carry a library suffix (<barcode>_<library>);
        # subset-bam matches the CB tag of the BAM, so the suffix is replaced by
        # the 10x "-1" suffix.
        tail -n +2 "${meta_file}" \
            | awk -v max="${DOUBLET_SCORE_MAX}" '$4 < max' \
            | cut -f1,5 \
            | grep "${donor}" \
            | cut -f1 \
            | sed "s/_${library}/-1/" > "${barcode_file}"

        # bamtofastq refuses to write into an existing directory
        if [ -d "${FASTQ_DIR}/${full_name}_FASTQ" ]; then
            rm -r "${FASTQ_DIR}/${full_name}_FASTQ"
        fi

        # --- Extract this donor's reads from the pooled BAM --------------------
        subset-bam \
            --bam "${CELLRANGER_DIR}/${group}/${library}/outs/possorted_genome_bam.bam" \
            --cell-barcodes "${barcode_file}" \
            --cores 16 \
            --out-bam "${FASTQ_DIR}/${full_name}.bam"

        # --- Convert back to FASTQ for re-alignment ---------------------------
        bamtofastq \
            --nthreads=16 \
            "${FASTQ_DIR}/${full_name}.bam" \
            "${FASTQ_DIR}/${full_name}_FASTQ"

    done < "${donor_list_dir}/${library}.txt"
done < "${META_DIR}/sample_list_group.tsv"

echo "Per-donor FASTQ files written to ${FASTQ_DIR}"
