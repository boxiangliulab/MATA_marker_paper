${cellranger} count --id=B1_L1 \
	--fastqs=./B1_L1 \
	--sample=B1_L1 \
	--transcriptome=./refdata-gex-GRCh38-2020-A \
	--include-introns true \
	--create-bam true
echo "finish ${id} cellranger count"
