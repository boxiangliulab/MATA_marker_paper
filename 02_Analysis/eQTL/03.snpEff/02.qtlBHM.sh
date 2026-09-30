#PBS -q normal
#PBS -l select=1:ncpus=16:mem=1000G
#PBS -l walltime=24:00:00
#PBS -P 11003054
#PBS -N qtlBHM
#PBS -o qtlBHM.o
#PBS -e qtlBHM.e

module load miniforge3
conda activate /data/projects/11003054/e1101943/biosoft/conda_envs/qtlbhm27
#cd /home/users/nus/e1101943/proj/qtlBHM

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(realpath "$SCRIPT_DIR/../../../../../../../..")"
QTLBHM_DIR="/data/projects/11003054/e1101943/biosoft/qtlBHM"  # set to your qtlBHM installation path
QTLBHM_PYTHON="/data/projects/11003054/e1101943/biosoft/conda_envs/qtlbhm27/bin/python"  # set to your qtlBHM conda python

${QTLBHM_PYTHON} ${QTLBHM_DIR}/infer_causal_variants.py \
  --output_prefix ${PROJECT_ROOT}/data/09.eQTL/03.Saige_QTL/02.output/05.snpEff/02.output/ALL.leadSNP.variant_effects.qtlBHM \
  ${PROJECT_ROOT}/data/09.eQTL/03.Saige_QTL/02.output/05.snpEff/02.output/input_for_qtlBHM_clean.tsv.gz \
  ${PROJECT_ROOT}/data/09.eQTL/03.Saige_QTL/02.output/05.snpEff/02.output/ALL.leadSNP.variant_effects.qtlBHM_input.bed.gz
