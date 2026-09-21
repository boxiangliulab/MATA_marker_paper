import os
import subprocess

traits = [
    "Acute_Pancreatitis_ldsc", "BMI_ldsc", "Body_Weight_ldsc", "Chronic_Pancreatitis_ldsc",
    "Coronary_Artery_Disease_ldsc", "Diabetic_Nephropathy_ldsc", "Glucocorticoids_ldsc",
    "Glucose_ldsc", "HbA1c_ldsc", "HDL_ldsc", "Height_ldsc", "LDL_ldsc",
    "Pancreatic_Cancer_ldsc", "Total_Cholesterol_ldsc", "Triglycerides_ldsc",
    "Waist_Circumference_ldsc", "2hGlu_ldsc", "Adiponectin_ldsc", "Fasting_Glucose_ldsc",
    "Fasting_Insulin_ldsc", "T2DSakaue_ldsc", "T2DSpracklen_ldsc",
]

sat_celltypes = [
    "Adipocyte", "Committed_preadipocyte", "CEC", "VEC", "AEC", "PVM",
    "SMC", "Early_preadipocyte", "LAM", "Areg", "Pericyte",
]

python_executable = "/home/users/nus/e1124850/anaconda3/envs/ldsc39/bin/python"
ldsc = "/home/users/nus/e1124850/scratch/github/ldsc/ldsc.py"
sumstats_dir = "/home/users/nus/e1124850/scratch"
annotation_dir = "/home/users/nus/e1124850/e1124850/co_lab/tianchi/gene_annot_file_grch37"
baseline_prefix = "/home/users/nus/e1124850/e1124850/co_lab/tongyihan/1000G_Phase3_EAS_baselineLD_v2.2_ldscores/baselineLD."
frequency_prefix = "/home/users/nus/e1124850/e1124850/co_lab/tongyihan/round1/1000G_Phase3_EAS_plinkfiles/1000G.EAS.QC."
weights_prefix = "/home/users/nus/e1124850/e1124850/co_lab/tongyihan/round1/1000G_Phase3_EAS_weights_hm3_no_MHC/weights.EAS.hm3_noMHC."
output_dir = "/home/users/nus/e1124850/scratch/sqtl_h2"

os.makedirs(output_dir, exist_ok=True)

run_env = os.environ.copy()
run_env["PYTHONNOUSERSITE"] = "1"
run_env.pop("PYTHONPATH", None)

failed_runs = []

for trait in traits:
    sumstats_file = f"{sumstats_dir}/CELLECT_{trait}.sumstats.gz"

    for celltype in sat_celltypes:
        annotation_prefix = f"{annotation_dir}/tianchi_sat_{celltype}_sqtlsig_gene_annot_file."
        output_prefix = f"{output_dir}/tianchi_sat_{trait}_{celltype}"

        cmd = [
            python_executable, "-s", ldsc,
            "--h2", sumstats_file,
            "--ref-ld-chr", f"{annotation_prefix},{baseline_prefix}",
            "--overlap-annot",
            "--frqfile-chr", frequency_prefix,
            "--w-ld-chr", weights_prefix,
            "--out", output_prefix,
            "--print-coefficients",
        ]

        print(f"\nRunning SAT: {trait} | {celltype}")

        try:
            subprocess.run(cmd, check=True, env=run_env)
            print(f"SUCCESS: {trait} | {celltype}")
        except subprocess.CalledProcessError as error:
            print(f"FAILED: {trait} | {celltype} | exit code={error.returncode}")
            failed_runs.append((trait, celltype, error.returncode))
        except Exception as error:
            print(f"UNEXPECTED ERROR: {trait} | {celltype} | {error}")
            failed_runs.append((trait, celltype, str(error)))

failure_file = f"{output_dir}/tianchi_sat_failed_runs.tsv"

with open(failure_file, "w") as handle:
    handle.write("trait\tcelltype\terror\n")
    for trait, celltype, error in failed_runs:
        handle.write(f"{trait}\t{celltype}\t{error}\n")

print(f"\nAll SAT jobs finished. Failed jobs: {len(failed_runs)}")
print(f"Failure list: {failure_file}")