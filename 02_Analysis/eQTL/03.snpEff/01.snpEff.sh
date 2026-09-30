# The reference used is same as CellRanger
# config your snpEff first: https://www.jianshu.com/p/a6e46d0c07ee
# path to snpEff: set SNPEFF_DIR below
module load miniforge3
conda activate seurat
#java -jar ~/snpEff/snpEff.jar build -c ~/snpEff/snpEff.config -gtf22 -v hg38

##### v1 ##########
#java -Xmx10G -jar ./snpEff.jar eff -c ./snpEff.config hg38 /home/users/nus/e1101943/scratch/project/00.adipose/data/09.eQTL/01.SAIGE_QTL/05.snpEff/SAT/unique_lead_eSNP.vcf.recode.vcf > /home/users/nus/e1101943/scratch/project/00.adipose/data/09.eQTL/01.SAIGE_QTL/05.snpEff/SAT/unique_lead_eSNP_snpEff_annotated.vcf -csvStats /home/users/nus/e1101943/scratch/project/00.adipose/data/09.eQTL/01.SAIGE_QTL/05.snpEff/SAT/unique_lead_eSNP_snpEff_annotated.csv -stats /home/users/nus/e1101943/scratch/project/00.adipose/data/09.eQTL/01.SAIGE_QTL/05.snpEff/SAT/unique_lead_eSNP_snpEff_annotated.html

#java -jar ./SnpSift.jar extractFields /home/users/nus/e1101943/scratch/project/00.adipose/data/09.eQTL/01.SAIGE_QTL/05.snpEff/SAT/unique_lead_eSNP_snpEff_annotated.vcf CHROM POS REF ALT "ANN[*].EFFECT" > /home/users/nus/e1101943/scratch/project/00.adipose/data/09.eQTL/01.SAIGE_QTL/05.snpEff/SAT/variant_effects_no_gene.tsv

##################
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(realpath "$SCRIPT_DIR/../../../../../../../..")"
SNPEFF_DIR="/data/projects/11003054/e1101943/biosoft/snpEff"  # set to your snpEff installation path

input="${PROJECT_ROOT}/data/09.eQTL/03.Saige_QTL/02.output/05.snpEff/01.input/"
output="${PROJECT_ROOT}/data/09.eQTL/03.Saige_QTL/02.output/05.snpEff/02.output"

java -Xmx10G -jar ${SNPEFF_DIR}/snpEff.jar eff -c ${SNPEFF_DIR}/snpEff.config hg38 $input/ALL.leadSNP.recode.vcf > $output/ALL.leadSNP.recode.snpEff_annotated.vcf -csvStats $output/ALL.leadSNP.snpEff_annotated.csv -stats $output/ALL.leadSNP.snpEff_annotated.html

java -jar ${SNPEFF_DIR}/SnpSift.jar extractFields -s "," -e "." $output/ALL.leadSNP.recode.snpEff_annotated.vcf CHROM POS REF ALT "ANN[*].EFFECT" > $output/ALL.leadSNP.variant_effects.tsv
