templete = '''BASE_OUTPUT_DIR: /home/users/nus/e1124850/e1124850/co_lab/tongyihan/round4/vat

SPECIFICITY_INPUT:
  - id: Adipose_VAT
    path: /home/users/nus/e1124850/e1124850/co_lab/tongyihan/round4/VAT_top_10percent.csv.gz

GWAS_SUMSTATS:
  - id: 2hGlu
    path: /home/users/nus/e1124850/scratch/CELLECT_2hGlu_ldsc.sumstats.gz
  - id: Acute_Pancreatitis
    path: /home/users/nus/e1124850/scratch/CELLECT_Acute_Pancreatitis_ldsc.sumstats.gz
  - id: Adiponectin
    path: /home/users/nus/e1124850/scratch/CELLECT_Adiponectin_ldsc.sumstats.gz
  - id: BMI
    path: /home/users/nus/e1124850/scratch/CELLECT_BMI_ldsc.sumstats.gz
  - id: Body_Weight
    path: /home/users/nus/e1124850/scratch/CELLECT_Body_Weight_ldsc.sumstats.gz
  - id: Chronic_Pancreatitis
    path: /home/users/nus/e1124850/scratch/CELLECT_Chronic_Pancreatitis_ldsc.sumstats.gz
  - id: Coronary_Artery_Disease
    path: /home/users/nus/e1124850/scratch/CELLECT_Coronary_Artery_Disease_ldsc.sumstats.gz
  - id: Diabetic_Nephropathy
    path: /home/users/nus/e1124850/scratch/CELLECT_Diabetic_Nephropathy_ldsc.sumstats.gz
  - id: Fasting_Glucose
    path: /home/users/nus/e1124850/scratch/CELLECT_Fasting_Glucose_ldsc.sumstats.gz
  - id: Fasting_Insulin
    path: /home/users/nus/e1124850/scratch/CELLECT_Fasting_Insulin_ldsc.sumstats.gz
  - id: Glucocorticoids
    path: /home/users/nus/e1124850/scratch/CELLECT_Glucocorticoids_ldsc.sumstats.gz
  - id: Glucose
    path: /home/users/nus/e1124850/scratch/CELLECT_Glucose_ldsc.sumstats.gz
  - id: HbA1c
    path: /home/users/nus/e1124850/scratch/CELLECT_HbA1c_ldsc.sumstats.gz
  - id: HDL
    path: /home/users/nus/e1124850/scratch/CELLECT_HDL_ldsc.sumstats.gz
  - id: Height
    path: /home/users/nus/e1124850/scratch/CELLECT_Height_ldsc.sumstats.gz
  - id: LDL
    path: /home/users/nus/e1124850/scratch/CELLECT_LDL_ldsc.sumstats.gz
  - id: Pancreatic_Cancer
    path: /home/users/nus/e1124850/scratch/CELLECT_Pancreatic_Cancer_ldsc.sumstats.gz
  - id: Total_Cholesterol
    path: /home/users/nus/e1124850/scratch/CELLECT_Total_Cholesterol_ldsc.sumstats.gz
  - id: Triglycerides
    path: /home/users/nus/e1124850/scratch/CELLECT_Triglycerides_ldsc.sumstats.gz
  - id: T2DSakaue
    path: /home/users/nus/e1124850/scratch/CELLECT_T2DSakaue_ldsc.sumstats.gz
  - id: T2DSpracklen
    path: /home/users/nus/e1124850/scratch/CELLECT_T2DSpracklen_ldsc.sumstats.gz
  - id: Waist_Circumference
    path: /home/users/nus/e1124850/scratch/CELLECT_Waist_Circumference_ldsc.sumstats.gz

ANALYSIS_TYPE: 
  prioritization: False 
  conditional: False 
  heritability: True 
  heritability_intervals: False 

WINDOW_DEFINITION: 
  WINDOW_SIZE_KB:
    100

GENE_COORD_FILE: 
  /home/users/nus/e1124850/CELLECT/data/shared/gene_coordinates.GRCh37.ensembl_v91.txt
KEEP_ANNOTS: 
  False


HERITABILITY_INPUT:
  - id: Adipose_VAT
    annotations: [ANNOTATIONS]

LDSC_CONST:
  DATA_DIR: 
    /home/users/nus/e1124850/CELLECT/data/ldsc/
  LDSC_DIR: 
    /home/users/nus/e1124850/ldsc/
  NUMPY_CORES:
    1

'''



celltypes = ['Adipocyte', 'AEC', 'Areg', 'B', 'CD4pos_T', 'CD56dim_CD16pos_NK', 'CD8pos_T', 'cDC2', 'CEC', 'Classical_monocyte', 'Committed_preadipocyte', 'Early_preadipocyte', 'IGFBP2pos_cell', 'LAM', 'LEC', 'Mesothelial', 'Pericyte', 'PVM', 'SMC', 'TIM4pos_ATM', 'TIM4pos_CD11cpos_ATM', 'VEC']

annotations = ', '.join("'{}'".format(c) for c in celltypes)

with open('/Users/theeeight/github/heritability_enrichment_ldsc/round4/vatconfig/config_vat.yaml', 'w') as f:
    f.write(templete.replace('ANNOTATIONS', annotations))
