#!/bin/bash
### Remove relateness
# Accept command-line argument
celltype=$1
depot=$2

# Project root (8 levels up from this script)
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(cd "${SCRIPT_DIR}/../../../../../../../.." && pwd)

# Define paths
# for data output
output="${PROJECT_ROOT}/data/09.eQTL/03.Saige_QTL/02.output/03.formal_calling/ancestry_specific_version/"
#output="/home/users/nus/e1101943/scratch/project/00.adipose/data/09.eQTL/01.SAIGE_QTL/03.formal_calling/"
mkdir -p "$output/$depot"
mkdir -p "$output/$depot/$celltype"
output="$output/$depot/$celltype"

# for script output
output_dir="${PROJECT_ROOT}/script/09.eQTL/03.Saige_QTL/04.QTL_calling/02.output/03.formal_calling/ancestry_specific_version/"
#output_dir="/home/users/nus/e1101943/scratch/project/00.adipose/script/01.adipose_sqtl/13.eQTL/03.Saige_QTL/02.output/03.formal_calling/"
mkdir -p "$output_dir/$depot"
mkdir -p "$output_dir/$depot/$celltype"
output_dir="$output_dir/$depot/$celltype"

batch_size=400
#phenotype="${PROJECT_ROOT}/data/09.eQTL/03.Saige_QTL/01.input/01.phenotype/${depot}/${celltype}/${celltype}_saige_phenotype.RNA.safe_gene_name.gPC8_pPC10.tsv" # SAT
phenotype="${PROJECT_ROOT}/data/09.eQTL/03.Saige_QTL/01.input/01.phenotype/${depot}/${celltype}/${celltype}_saige_phenotype.RNA.safe_gene_name.gPC6_pPC10.tsv" # VAT
step1_genotype="${PROJECT_ROOT}/data/09.eQTL/03.Saige_QTL/01.input/03.step01_plink/ancestry_specific/adipose_after_imputation_296_plink.pruned.mac20.random2000"

# Check if input file exists
if [ ! -f "$phenotype" ]; then
    echo "Input file does not exist: $phenotype"
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Initialize batch counters
counter=0
batch_num=1

# Read genes from the 12th column onward and generate scripts
#awk -F'\t' 'NR==1 {for (i=26; i<=NF; i++) print $i}' "$phenotype" | while IFS= read -r gene; do #SAT
awk -F'\t' 'NR==1 {for (i=24; i<=NF; i++) print $i}' "$phenotype" | while IFS= read -r gene; do #VAT

    # Skip empty lines
    if [[ -z "$gene" ]]; then
        continue
    fi

    # Start a new batch script if counter is 0
    if (( counter == 0 )); then
        script_file="${output_dir}/batch${batch_num}.sh"
        echo "#!/bin/bash" > "$script_file"
        echo "#PBS -q normal" >> "$script_file"
        echo "#PBS -l select=1:ncpus=2:mem=100G" >> "$script_file"
        echo "#PBS -l walltime=24:00:00" >> "$script_file"
        echo "#PBS -P 11003054" >> "$script_file"
        echo "#PBS -N Saige_step01_${celltype}_batch${batch_num}" >> "$script_file"
        echo "#PBS -o Saige_step01_${celltype}_batch${batch_num}.o" >> "$script_file"
        echo "#PBS -e Saige_step01_${celltype}_batch${batch_num}.e" >> "$script_file"
        echo "module load miniforge3" >> "$script_file"
        echo "conda activate sr" >> "$script_file"
        cd ${PBS_O_WORKDIR}
        echo "" >> "$script_file"
    fi

    # step01
    echo "echo \"[step01] Processing gene: $gene\"" >> "$script_file"
    echo "/data/projects/11003054/e1101943/biosoft/conda_envs/saigeqtl/bin/Rscript ${PROJECT_ROOT}/script/09.eQTL/01.AIDA_eQTL/01.SAIGE_QTL/step1_fitNULLGLMM_qtl.R \\" >> "$script_file"
    echo "--useSparseGRMtoFitNULL=FALSE \\" >> "$script_file"
    echo "--useGRMtoFitNULL=FALSE \\" >> "$script_file"
    echo "--phenoFile=$phenotype \\" >> "$script_file"
    echo "--phenoCol=$gene \\" >> "$script_file"
    #echo "--covarColList=sex,diabetic,Age,BMI,log1p_UMI,percent.mt,gPC1,gPC2,gPC3,gPC4,gPC5,gPC6,gPC7,gPC8,pPC1,pPC2,pPC3,pPC4,pPC5,pPC6,pPC7,pPC8,pPC9,pPC10 \\" >> "$script_file" # SAT
    #echo "--sampleCovarColList=sex,diabetic,Age,BMI,gPC1,gPC2,gPC3,gPC4,gPC5,gPC6,gPC7,gPC8 \\" >> "$script_file" # SAT
    echo "--covarColList=sex,diabetic,Age,BMI,log1p_UMI,percent.mt,gPC1,gPC2,gPC3,gPC4,gPC5,gPC6,pPC1,pPC2,pPC3,pPC4,pPC5,pPC6,pPC7,pPC8,pPC9,pPC10 \\" >> "$script_file" # VAT
    echo "--sampleCovarColList=sex,diabetic,Age,BMI,gPC1,gPC2,gPC3,gPC4,gPC5,gPC6 \\" >> "$script_file" #VAT
    echo "--sampleIDColinphenoFile=donor_id \\" >> "$script_file"
    echo "--traitType=count \\" >> "$script_file"
    echo "--outputPrefix=$output/$gene \\" >> "$script_file"
    echo "--plinkFile=$step1_genotype \\" >> "$script_file"
    echo "--skipVarianceRatioEstimation=FALSE \\" >> "$script_file"
    echo "--isRemoveZerosinPheno=FALSE \\" >> "$script_file"
    echo "--isCovariateOffset=FALSE \\" >> "$script_file"
    echo "--tol=0.00001 \\" >> "$script_file"
    echo "--skipModelFitting=FALSE \\"  >> "$script_file"
    echo "--maxiterPCG=5000 \\"  >> "$script_file"
    echo "--nThreads=16 \\"  >> "$script_file"
    echo "--IsOverwriteVarianceRatioFile=TRUE" >> "$script_file"
    echo "" >> "$script_file"
    
    # Increment counters
    ((counter++))

    # Check if batch size is reached
    if (( counter == batch_size )); then
        counter=0
        ((batch_num++))
    fi
done

# Handle leftover genes for the last batch
if (( counter > 0 )); then
    echo "Final batch contains $counter genes. Scripts have been generated."
fi

echo "Script generation complete. All scripts are saved in: $output_dir"

