#PBS -q normal
#PBS -l select=1:ncpus=1:mem=300G
#PBS -l walltime=2:00:00
#PBS -N imputation_results
#PBS -o imputation_results.o
#PBS -e imputation_results.e
#PBS -P 11003054

module load parallel

WORKDIR="/project/adipose/fullgenotype/new_version/14.imputation_results"
cd $WORKDIR

curl -sL https://imputationserver.sph.umich.edu/get/wJJOw9eaTaHDAyF0w78yTQzFRhaINOKHkIClQLUb | bash

ls *.zip | parallel -j+0 '~/miniconda3/envs/qtl/bin/7z x -p"VWQF+Nzd6a9jw}" {}'
