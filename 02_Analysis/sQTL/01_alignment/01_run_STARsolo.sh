#!/bin/bash
# Per-donor STARsolo alignment with WASP filtering.
#
# Re-align the demultiplexed, donor-specific FASTQ files with STARsolo. Each
# donor is aligned against its own genotypes (--waspOutputMode SAMtag with
#   --varVCFfile) so that reads showing allelic mapping bias can be removed
# before junction quantification, which is essential for unbiased sQTL
# mapping.
#
# Library chemistry: 10x 5' gene expression v3
#   - 16 bp cell barcode + 12 bp UMI on read 1 (--soloBarcodeMate 1)
#   - the first 41 bases of read 1 (barcode + UMI + TSO) are clipped

set -euo pipefail
cd "$(dirname "$0")"

# Paths
META_DIR=../data/00_metadata
FASTQ_DIR=../data/00_bamtofastq
STARSOLO_DIR=../data/01_STARsolo
# Per-donor VCF used for WASP filtering
GENOTYPE_BY_DONOR_DIR=../data/genotype/vcf_by_individual
# Reference: STAR index, annotation, and the 10x 5' GEX v3 barcode whitelist
STAR_INDEX=../data/reference/STAR_index/hg38_gencode_v32
GTF_FILE=../data/reference/gencode.v32.primary_assembly.annotation.gtf
CB_WHITELIST=../data/reference/3M-5pgex-jan-2023.txt

mkdir -p "${STARSOLO_DIR}"

while IFS=$'\t' read -r group library; do
    echo "=== ${group} / ${library} ==="

    while read -r donor; do
        echo "--- donor ${donor} ---"
        full_name="${library}_${donor}"
        out_dir="${STARSOLO_DIR}/${full_name}"
        mkdir -p "${out_dir}/log_files"

        # bamtofastq splits the reads over several chunks; STAR accepts them as
        # a comma-separated list.
        R1=$(ls "${FASTQ_DIR}/${full_name}_FASTQ"/*/*R1* | tr '\n' ',' | sed 's/,$//')
        R2=$(ls "${FASTQ_DIR}/${full_name}_FASTQ"/*/*R2* | tr '\n' ',' | sed 's/,$//')

        STAR \
            --runThreadN 24 \
            --limitBAMsortRAM 1024000000000 \
            --genomeDir "${STAR_INDEX}" \
            --twopassMode Basic \
            --soloStrand Forward \
            --soloType CB_UMI_Simple \
            --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 \
            --soloBarcodeMate 1 --clip5pNbases 41 0 \
            --soloCBwhitelist "${CB_WHITELIST}" \
            --soloCBmatchWLtype 1MM \
            --soloUMIdedup 1MM_Directional_UMItools \
            --waspOutputMode SAMtag --varVCFfile "${GENOTYPE_BY_DONOR_DIR}/${donor}.vcf" \
            --readFilesIn "${R1}" "${R2}" \
            --readFilesCommand zcat \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM XS vG vA vW \
            --sjdbGTFfile "${GTF_FILE}" \
            --soloOutFileNames "${full_name}" \
            --limitSjdbInsertNsj 2000000 \
            --outFileNamePrefix "${out_dir}/${full_name}_" \
            1> "${out_dir}/log_files/${donor}.out" \
            2> "${out_dir}/log_files/${donor}.err"

    done < "${META_DIR}/donor_list/${library}.txt"
done < "${META_DIR}/sample_list_group.tsv"

# Manifest of aligned samples (one "<library>_<donor>" per line), used by the
# junction-extraction stage.
ls -d "${STARSOLO_DIR}"/*/ | xargs -n 1 basename | sort -V > "${META_DIR}/list_samples.txt"
echo "Alignment finished; sample list written to ${META_DIR}/list_samples.txt"
