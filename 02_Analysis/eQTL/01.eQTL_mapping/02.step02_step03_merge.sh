#!/bin/bash
### Remove relateness
#### Initially, save the file in the project, but disk quota exceed, so I change the save path to scratch
# Accept command-line argument
celltype=$1
depot=$2

# Project root (8 levels up from this script)
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(cd "${SCRIPT_DIR}/../../../../../../../.." && pwd)

# Define paths
output="${PROJECT_ROOT}/data/09.eQTL/03.Saige_QTL/02.output/03.formal_calling/ancestry_specific_version/"
#output="/home/users/nus/e1101943/scratch/project/00.adipose/script/01.adipose_sqtl/13.eQTL/03.Saige_QTL/02.output/03.formal_calling/ancestry_specific_version/"
mkdir -p "$output/$depot"
mkdir -p "$output/$depot/$celltype"
output="$output/$depot/$celltype"

output_dir="${PROJECT_ROOT}/script/09.eQTL/03.Saige_QTL/04.QTL_calling/02.output/03.formal_calling/ancestry_specific_version/"
#output_dir="/home/users/nus/e1101943/scratch/project/00.adipose/script/01.adipose_sqtl/13.eQTL/03.Saige_QTL/02.output/03.formal_calling"
mkdir -p "$output_dir/$depot"
mkdir -p "$output_dir/$depot/$celltype"
mkdir -p "$output_dir/$depot/$celltype/step02_03"
output_dir="$output_dir/$depot/$celltype/step02_03"

batch_size=400
phenotype="${PROJECT_ROOT}/data/09.eQTL/03.Saige_QTL/01.input/01.phenotype/$depot/$celltype/${celltype}_saige_phenotype.RNA.safe_gene_name.gPC8_pPC10.tsv" # SAT
#phenotype="${PROJECT_ROOT}/data/09.eQTL/03.Saige_QTL/01.input/01.phenotype/$depot/$celltype/${celltype}_saige_phenotype.RNA.safe_gene_name.gPC6_pPC10.tsv" # VAT
bed_file="${PROJECT_ROOT}/data/09.eQTL/02.saige_vs_QTLtools/gene_positions_safe_name.bed"

# Build hash pair of gene and chr
declare -A gene_chr
while IFS=$'\t' read -r chr _ _ _ _ gene _; do
    gene_chr["$gene"]=$chr
done < "$bed_file"

# Check if input file exists
if [ ! -f "$phenotype" ]; then
    echo "Input file does not exist: $phenotype" >&2
    exit 1
fi

# Initialize batch counters
counter=0
batch_num=1

# Process genes
while IFS= read -r gene; do
    # Skip empty lines
    [[ -z "$gene" ]] && continue
    
    # Check if gene exists in BED file
    if [[ -z "${gene_chr[$gene]}" ]]; then
        echo "WARNING: Gene $gene not found in BED file" >&2
        continue
    fi
    
    chr="${gene_chr[$gene]}"
    echo "Processing gene: $gene (Chromosome: $chr)"
    
    # Start a new batch script if counter is 0
    if (( counter == 0 )); then
        script_file="${output_dir}/step02_03_batch${batch_num}.sh"
        cat > "$script_file" <<EOF
#!/bin/bash
#PBS -q normal
#PBS -l select=1:ncpus=4:mem=300G
#PBS -l walltime=24:00:00
#PBS -P 11003054
#PBS -N Saige_step02_${celltype}_batch${batch_num}
#PBS -o ${output_dir}/Saige_step02_${celltype}_batch${batch_num}.o
#PBS -e ${output_dir}/Saige_step02_${celltype}_batch${batch_num}.e

cd \${PBS_O_WORKDIR}

EOF
    fi
    
    # step02 commands
    cat >> "$script_file" <<EOF
echo "[step02] Processing gene: $gene ($chr)"
/data/projects/11003054/e1101943/biosoft/conda_envs/saigeqtl/bin/Rscript \\
    ${PROJECT_ROOT}/script/09.eQTL/01.AIDA_eQTL/01.SAIGE_QTL/step2_tests_qtl.R \\
    --vcfFile="/home/project/11003054/e1353486/project/adipose/fullgenotype/Ancestry_specific/02.vcf_by_chr_edited_rsid/adipose_after_imputation_296_merged_genome_rename_rsID_nochr_neur_${chr}.vcf.gz" \\
    --vcfFileIndex="/home/project/11003054/e1353486/project/adipose/fullgenotype/Ancestry_specific/02.vcf_by_chr_edited_rsid/adipose_after_imputation_296_merged_genome_rename_rsID_nochr_neur_${chr}.vcf.gz.csi" \\
    --vcfField="GT" \\
    --SAIGEOutputFile="$output/${gene}_cis" \\
    --chrom="${chr}" \\
    --minMAF=0.05 \\
    --minMAC=20 \\
    --LOCO=FALSE \\
    --GMMATmodelFile="$output/$gene.rda" \\
    --SPAcutoff=2 \\
    --varianceRatioFile="$output/$gene.varianceRatio.txt" \\
    --markers_per_chunk=20000 \\
    --rangestoIncludeFile="${PROJECT_ROOT}/data/09.eQTL/02.saige_vs_QTLtools/cis_region_file_for_each_gene/$gene.tsv"

EOF
    
    # step03 commands
    cat >> "$script_file" <<EOF
echo "[step03] Processing gene: $gene ($chr)"
/data/projects/11003054/e1101943/biosoft/conda_envs/saigeqtl/bin/Rscript \\
    ${PROJECT_ROOT}/script/09.eQTL/01.AIDA_eQTL/01.SAIGE_QTL/step3_gene_pvalue_qtl.R \\
    --assocFile="$output/${gene}_cis" \\
    --geneName="$gene" \\
    --genePval_outputFile="$output/${gene}_cis_genePval"

EOF

    ((counter++))
    
    # Check if batch size is reached
    if (( counter == batch_size )); then
        counter=0
        ((batch_num++))
    fi
done < <(awk -F'\t' 'NR==1 {for (i=26; i<=NF; i++) print $i}' "$phenotype") # modify this, VAT = 24 and SAT starts from 26

# Handle leftover genes for the last batch
if (( counter > 0 )); then
    echo "Final batch contains $counter genes. Scripts have been generated."
fi

echo "Script generation complete. All scripts are saved in: $output_dir"
