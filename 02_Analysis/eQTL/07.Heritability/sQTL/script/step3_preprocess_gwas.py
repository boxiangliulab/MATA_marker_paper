

import pandas as pd
import os

gwas_home_path = '/data/projects/11003054/e1101943/project/01.adipose_sqtl/data/27.GWAS/prepare_required_columns/with_rsid'

gwas_file_name_dict = {
    # '2hGlu':'2hGlu_Chen_2021_Nature_Genetics_hg38.txt', ## 
    'Acute_Pancreatitis': 'Acute_Pancreatitis_Sakaue_2021_Nature_Genetics_hg38.txt',
    # 'Adiponectin': 'Adiponectin_Wu_2014_Human_Molecular_Genetics_hg38.txt', ## 
    'BMI': 'BMI_Sakaue_2021_Nature_Genetics_hg38.txt',
    "Body_Weight":'Body_Weight_Sakaue_2021_Nature_Genetics_hg38.txt',
    'Chronic_Pancreatitis':'Chronic_Pancreatitis_Sakaue_2021_Nature_Genetics_hg38.txt',
    'Coronary_Artery_Disease': 'Coronary_Artery_Disease_Sakaue_2021_Nature_Genetics_hg38.txt',
    'Diabetic_Nephropathy': 'Diabetic_Nephropathy_Sakaue_2021_Nature_Genetics_hg38.txt',
    # 'Fasting_Glucose': 'Fasting_Glucose_Chen_2021_Nature_Genetics_hg38.txt', ## 
    # 'Fasting_Insulin': 'Fasting_Insulin_Chen_2021_Nature_Genetics_hg38.txt', ## 
    'Glucocorticoids': 'Glucocorticoids_Sakaue_2021_Nature_Genetics_hg38.txt',
    'Glucose': 'Glucose_Sakaue_2021_Nature_Genetics_hg38.txt',
    'HbA1c': 'HbA1c_Sakaue_2021_Nature_Genetics_hg38.txt',
    'HDL': 'HDL_Cholesterol_Sakaue_2021_Nature_Genetics_hg38.txt',
    "Height":'Height_Sakaue_2021_Nature_Genetics_hg38.txt',
    'LDL': 'LDL_Cholesterol_Sakaue_2021_Nature_Genetics_hg38.txt',
    'Pancreatic_Cancer': 'Pancreatic_Cancer_Sakaue_2021_Nature_Genetics_hg38.txt',
    'Total_Cholesterol': 'Total_Cholesterol_Sakaue_2021_Nature_Genetics_hg38.txt',
    'Triglycerides': 'Triglycerides_Sakaue_2021_Nature_Genetics_hg38.txt',
    # 'T2DSakaue': 'Type2_Diabetes_Sakaue_2021_Nature_Genetics_hg38.txt', ## 
    # 'T2DSpracklen': 'Type2_Diabetes_Spracklen_2020_Nature_hg38.txt', ##
    "Waist_Circumference":'Waist_Circumference_Nam_2022_Cell_Genomics_hg38.txt'
    }

gwas_samplesize_dict = {
    '2hGlu': 8509,
    'Acute_Pancreatitis': 178298,
    'Adiponectin': 7825,
    'BMI': 163835,
    "Body_Weight": 165419,
    'Chronic_Pancreatitis': 177928,
    'Coronary_Artery_Disease': 178726,
    'Diabetic_Nephropathy': 132984,
    'Fasting_Glucose': 35619,
    'Fasting_Insulin': 29792,
    'Glucocorticoids': 178726,
    'Glucose': 133336,
    'HbA1c': 71221,
    'HDL': 74970,
    "Height": 165056,
    'LDL': 72866,
    'Pancreatic_Cancer': 159700,
    'Total_Cholesterol': 135808,
    'Triglycerides': 111667,
    'T2DSakaue': 177415,
    'T2DSpracklen': 433540,
    "Waist_Circumference": 72222}

for gwas in gwas_file_name_dict.keys():
    gwas_file_path = os.path.join(gwas_home_path, gwas_file_name_dict[gwas])
    gwas_df = pd.read_csv(gwas_file_path, sep='\t')
    print(f'{gwas} columns: {gwas_df.columns}')
    if len(gwas_df.columns) == 1:
        gwas_df = pd.read_csv(gwas_file_path, sep=' ')
    gwas_df = gwas_df.dropna()
    gwas_df['N'] = gwas_samplesize_dict[gwas]
    gwas_df['A1'] = gwas_df['SNPID_hg19'].apply(lambda x: x.split(':')[-1]) # effect allele
    gwas_df['A2'] = gwas_df['SNPID_hg19'].apply(lambda x: x.split(':')[-2]) # non-effect allele
    gwas_df['SNP'] = 'chr' + gwas_df['SNPID_hg38'] + ':' + gwas_df['A2'] + ':' + gwas_df['A1']
    if 'rs_id' in gwas_df.columns:
        gwas_df['SNP'] = gwas_df['rs_id']
    elif 'rsID' in gwas_df.columns:
        gwas_df['SNP'] = gwas_df['rsID']
    # gwas_df['SNP'] = gwas_df['rs_id'] / rsID
    gwas_df['Z'] = gwas_df['BETA'] / gwas_df['SE']
    gwas_df = gwas_df[['SNP', 'N', 'Z', 'A1', 'A2','Pvalue']]
    # gwas_df = gwas_df[['SNP', 'N', 'BETA', 'A1', 'A2','Pvalue']] for "Diabetic_Nephropathy"
    gwas_df.columns = ['SNP', 'N', 'Z', 'A1', 'A2','P']
    output_path = f'/home/users/nus/e1124850/scratch/{gwas}.ldsc.format.tsv'
    gwas_df.to_csv(output_path, sep='\t', index=False)


for gwas in gwas_file_name_dict.keys():
    output_path = f'/home/users/nus/e1124850/scratch/{gwas}.ldsc.format.tsv'
    print(f"python /home/users/nus/e1124850/ldsc/munge_sumstats.py --sumstats {output_path} --out /home/users/nus/e1124850/scratch/{gwas}_ldsc")


# python munge_sumstats.py --sumstats /home/users/nus/e1124850/scratch/Waist_Circumference.ldsc.format.tsv --out Waist_CircumferenceWaist_Circumference


'''
#PBS -q normal
#PBS -l select=1:ncpus=2:mem=300G
#PBS -l walltime=24:00:00
#PBS -P 11003054
#PBS -N step3_preprocess_gwas
#PBS -o /home/users/nus/e1124850/scratch/qsub_dir/step3_preprocess_gwas.o
#PBS -e /home/users/nus/e1124850/scratch/qsub_dir/step3_preprocess_gwas.e

source /home/project/11003054/e1124850/MM/software/anaconda3/etc/profile.d/conda.sh
conda deactivate
conda deactivate
conda activate bio

python /home/users/nus/e1124850/scratch/qsub_dir/step3_preprocess_gwas.py


'''




'''
#PBS -q normal
#PBS -l select=1:ncpus=2:mem=300G
#PBS -l walltime=24:00:00
#PBS -P 11003054
#PBS -N munge_sumstats
#PBS -o /home/users/nus/e1124850/scratch/qsub_dir/munge_sumstats.o
#PBS -e /home/users/nus/e1124850/scratch/qsub_dir/munge_sumstats.e

source /home/project/11003054/e1124850/MM/software/anaconda3/etc/profile.d/conda.sh
conda deactivate
conda deactivate
conda activate ldsc39

# python /home/users/nus/e1124850/ldsc/munge_sumstats.py --sumstats /home/users/nus/e1124850/scratch/2hGlu.ldsc.format.tsv --out /home/users/nus/e1124850/scratch/2hGlu_ldsc
# python /home/users/nus/e1124850/ldsc/munge_sumstats.py --sumstats /home/users/nus/e1124850/scratch/Adiponectin.ldsc.format.tsv --out /home/users/nus/e1124850/scratch/Adiponectin_ldsc
# python /home/users/nus/e1124850/ldsc/munge_sumstats.py --sumstats /home/users/nus/e1124850/scratch/Fasting_Glucose.ldsc.format.tsv --out /home/users/nus/e1124850/scratch/Fasting_Glucose_ldsc
# python /home/users/nus/e1124850/ldsc/munge_sumstats.py --sumstats /home/users/nus/e1124850/scratch/Fasting_Insulin.ldsc.format.tsv --out /home/users/nus/e1124850/scratch/Fasting_Insulin_ldsc
# python /home/users/nus/e1124850/ldsc/munge_sumstats.py --sumstats /home/users/nus/e1124850/scratch/T2DSakaue.ldsc.format.tsv --out /home/users/nus/e1124850/scratch/T2DSakaue_ldsc
# python /home/users/nus/e1124850/ldsc/munge_sumstats.py --sumstats /home/users/nus/e1124850/scratch/T2DSpracklen.ldsc.format.tsv --out /home/users/nus/e1124850/scratch/T2DSpracklen_ldsc

python /home/users/nus/e1124850/ldsc/munge_sumstats.py --sumstats /home/users/nus/e1124850/scratch/Acute_Pancreatitis.ldsc.format.tsv --out /home/users/nus/e1124850/scratch/Acute_Pancreatitis_ldsc &
python /home/users/nus/e1124850/ldsc/munge_sumstats.py --sumstats /home/users/nus/e1124850/scratch/BMI.ldsc.format.tsv --out /home/users/nus/e1124850/scratch/BMI_ldsc &
python /home/users/nus/e1124850/ldsc/munge_sumstats.py --sumstats /home/users/nus/e1124850/scratch/Body_Weight.ldsc.format.tsv --out /home/users/nus/e1124850/scratch/Body_Weight_ldsc &
python /home/users/nus/e1124850/ldsc/munge_sumstats.py --sumstats /home/users/nus/e1124850/scratch/Chronic_Pancreatitis.ldsc.format.tsv --out /home/users/nus/e1124850/scratch/Chronic_Pancreatitis_ldsc &
python /home/users/nus/e1124850/ldsc/munge_sumstats.py --sumstats /home/users/nus/e1124850/scratch/Coronary_Artery_Disease.ldsc.format.tsv --out /home/users/nus/e1124850/scratch/Coronary_Artery_Disease_ldsc  &
wait


python /home/users/nus/e1124850/ldsc/munge_sumstats.py --sumstats /home/users/nus/e1124850/scratch/Diabetic_Nephropathy.ldsc.format.tsv --out /home/users/nus/e1124850/scratch/Diabetic_Nephropathy_ldsc &
python /home/users/nus/e1124850/ldsc/munge_sumstats.py --sumstats /home/users/nus/e1124850/scratch/Glucocorticoids.ldsc.format.tsv --out /home/users/nus/e1124850/scratch/Glucocorticoids_ldsc &
python /home/users/nus/e1124850/ldsc/munge_sumstats.py --sumstats /home/users/nus/e1124850/scratch/Glucose.ldsc.format.tsv --out /home/users/nus/e1124850/scratch/Glucose_ldsc &
python /home/users/nus/e1124850/ldsc/munge_sumstats.py --sumstats /home/users/nus/e1124850/scratch/HbA1c.ldsc.format.tsv --out /home/users/nus/e1124850/scratch/HbA1c_ldsc &
python /home/users/nus/e1124850/ldsc/munge_sumstats.py --sumstats /home/users/nus/e1124850/scratch/HDL.ldsc.format.tsv --out /home/users/nus/e1124850/scratch/HDL_ldsc &
python /home/users/nus/e1124850/ldsc/munge_sumstats.py --sumstats /home/users/nus/e1124850/scratch/Height.ldsc.format.tsv --out /home/users/nus/e1124850/scratch/Height_ldsc &
wait



python /home/users/nus/e1124850/ldsc/munge_sumstats.py --sumstats /home/users/nus/e1124850/scratch/LDL.ldsc.format.tsv --out /home/users/nus/e1124850/scratch/LDL_ldsc &
python /home/users/nus/e1124850/ldsc/munge_sumstats.py --sumstats /home/users/nus/e1124850/scratch/Pancreatic_Cancer.ldsc.format.tsv --out /home/users/nus/e1124850/scratch/Pancreatic_Cancer_ldsc &
python /home/users/nus/e1124850/ldsc/munge_sumstats.py --sumstats /home/users/nus/e1124850/scratch/Total_Cholesterol.ldsc.format.tsv --out /home/users/nus/e1124850/scratch/Total_Cholesterol_ldsc &
python /home/users/nus/e1124850/ldsc/munge_sumstats.py --sumstats /home/users/nus/e1124850/scratch/Triglycerides.ldsc.format.tsv --out /home/users/nus/e1124850/scratch/Triglycerides_ldsc &
python /home/users/nus/e1124850/ldsc/munge_sumstats.py --sumstats /home/users/nus/e1124850/scratch/Waist_Circumference.ldsc.format.tsv --out /home/users/nus/e1124850/scratch/Waist_Circumference_ldsc & 
wait



'''





# python /home/users/nus/e1124850/ldsc/munge_sumstats.py --sumstats /home/users/nus/e1124850/scratch/Acute_Pancreatitis.ldsc.format.tsv --out /home/users/nus/e1124850/scratch/Acute_Pancreatitis_ldsc
# python /home/users/nus/e1124850/ldsc/munge_sumstats.py --sumstats /home/users/nus/e1124850/scratch/BMI.ldsc.format.tsv --out /home/users/nus/e1124850/scratch/BMI_ldsc
# python /home/users/nus/e1124850/ldsc/munge_sumstats.py --sumstats /home/users/nus/e1124850/scratch/Body_Weight.ldsc.format.tsv --out /home/users/nus/e1124850/scratch/Body_Weight_ldsc
# python /home/users/nus/e1124850/ldsc/munge_sumstats.py --sumstats /home/users/nus/e1124850/scratch/Chronic_Pancreatitis.ldsc.format.tsv --out /home/users/nus/e1124850/scratch/Chronic_Pancreatitis_ldsc
# python /home/users/nus/e1124850/ldsc/munge_sumstats.py --sumstats /home/users/nus/e1124850/scratch/Coronary_Artery_Disease.ldsc.format.tsv --out /home/users/nus/e1124850/scratch/Coronary_Artery_Disease_ldsc
# python /home/users/nus/e1124850/ldsc/munge_sumstats.py --sumstats /home/users/nus/e1124850/scratch/Diabetic_Nephropathy.ldsc.format.tsv --out /home/users/nus/e1124850/scratch/Diabetic_Nephropathy_ldsc
# python /home/users/nus/e1124850/ldsc/munge_sumstats.py --sumstats /home/users/nus/e1124850/scratch/Glucocorticoids.ldsc.format.tsv --out /home/users/nus/e1124850/scratch/Glucocorticoids_ldsc
# python /home/users/nus/e1124850/ldsc/munge_sumstats.py --sumstats /home/users/nus/e1124850/scratch/Glucose.ldsc.format.tsv --out /home/users/nus/e1124850/scratch/Glucose_ldsc
# python /home/users/nus/e1124850/ldsc/munge_sumstats.py --sumstats /home/users/nus/e1124850/scratch/HbA1c.ldsc.format.tsv --out /home/users/nus/e1124850/scratch/HbA1c_ldsc
# python /home/users/nus/e1124850/ldsc/munge_sumstats.py --sumstats /home/users/nus/e1124850/scratch/HDL.ldsc.format.tsv --out /home/users/nus/e1124850/scratch/HDL_ldsc
# python /home/users/nus/e1124850/ldsc/munge_sumstats.py --sumstats /home/users/nus/e1124850/scratch/Height.ldsc.format.tsv --out /home/users/nus/e1124850/scratch/Height_ldsc
# python /home/users/nus/e1124850/ldsc/munge_sumstats.py --sumstats /home/users/nus/e1124850/scratch/LDL.ldsc.format.tsv --out /home/users/nus/e1124850/scratch/LDL_ldsc
# python /home/users/nus/e1124850/ldsc/munge_sumstats.py --sumstats /home/users/nus/e1124850/scratch/Pancreatic_Cancer.ldsc.format.tsv --out /home/users/nus/e1124850/scratch/Pancreatic_Cancer_ldsc
# python /home/users/nus/e1124850/ldsc/munge_sumstats.py --sumstats /home/users/nus/e1124850/scratch/Total_Cholesterol.ldsc.format.tsv --out /home/users/nus/e1124850/scratch/Total_Cholesterol_ldsc
# python /home/users/nus/e1124850/ldsc/munge_sumstats.py --sumstats /home/users/nus/e1124850/scratch/Triglycerides.ldsc.format.tsv --out /home/users/nus/e1124850/scratch/Triglycerides_ldsc
# python /home/users/nus/e1124850/ldsc/munge_sumstats.py --sumstats /home/users/nus/e1124850/scratch/Waist_Circumference.ldsc.format.tsv --out /home/users/nus/e1124850/scratch/Waist_Circumference_ldsc
