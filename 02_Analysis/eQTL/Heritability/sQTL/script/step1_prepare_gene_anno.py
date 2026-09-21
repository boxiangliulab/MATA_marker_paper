import pandas as pd
import os 
from glob import glob

from importlib import metadata
import pyranges as pr

## gene-coord-file
def __prepare_genecode(gencode):
    genecode_df = pr.read_gtf(gencode, as_df=True)
    genecode_df.drop(labels=genecode_df[genecode_df['Feature'] != 'gene'].index, inplace=True)
    genecode_df.reset_index(drop=True, inplace=True)
    genecode_df.drop(columns=[col for col in genecode_df if
                                col not in ['Chromosome', 'Start', 'End', 'Strand', 'gene_id', 'gene_name']],
                        inplace=True)
    genecode_df.rename({'Chromosome': 'chr', 'Start': 'start', 'End': 'end', 'Strand': 'strand'},
                        axis='columns', inplace=True)
    # gene_id_df = genecode_df['pheno_id'].str.split('.', n=1, expand=True)
    # genecode_df.loc[:, 'pheno_id'] = gene_id_df[0]

    genecode_df.loc[:, 'position'] = pd.NA
    genecode_df['position'].mask(genecode_df['strand'] == '+', genecode_df['start'], inplace=True)
    genecode_df['position'].mask(genecode_df['position'].isna(), genecode_df['end'], inplace=True)
    # genecode_df.drop(columns=['start', 'end'], inplace=True)
    genecode_df['chr'] = genecode_df['chr'].apply(lambda x: x.strip('chr'))
    return genecode_df


# gtf = '/data/projects/11003054/e1101943/project/00.reference/hg38/cellranger_reference/refdata-gex-GRCh38-2020-A/genes/genes.gtf'
# vcf = '/home/project/11003054/e1353486/project/adipose/fullgenotype/Ancestry_specific/00.merge_genotype/adipose_after_imputation_296_merged_genome_rename_rsID_nochr_neur.vcf.gz'

# gene_coord_file_fullcolumn = __prepare_genecode(gtf)
# gene_coord_file = gene_coord_file_fullcolumn[['gene_name','chr','start','end']]
# gene_coord_file.columns = ['GENE', 'CHR', 'START', 'END']
# gene_coord_file_path = '/home/users/nus/e1124850/scratch/github/heritability_enrichment_ldsc/round5sqtl/gene_coord_file.txt'
# gene_coord_file.to_csv(gene_coord_file_path, sep='\t', index=False)


## gene-set-file

sqtl_file = '/home/project/11003054/e1101920/adipose/data_v0/05_model/formal_analysis/results_new/All_sQTL_SAT_and_VAT_with_beta.tsv'
sqtl_df = pd.read_csv(sqtl_file, sep='\t')
sqtl_df = sqtl_df.dropna()
sat_sqtl_df = sqtl_df[sqtl_df['depot'] == 'SAT']
vat_sqtl_df = sqtl_df[sqtl_df['depot'] == 'VAT']


for celltype in sat_sqtl_df['celltype'].unique():
    celltype_df = sat_sqtl_df[sat_sqtl_df['celltype'] == celltype]
    celltype_df['gene'].drop_duplicates().to_csv(f'/home/users/nus/e1124850/scratch/github/heritability_enrichment_ldsc/round5sqtl/tianchi_sat_{celltype}_gene_set_file.txt', index=False, header=False)


for celltype in vat_sqtl_df['celltype'].unique():
    celltype_df = vat_sqtl_df[vat_sqtl_df['celltype'] == celltype]
    celltype_df['gene'].drop_duplicates().to_csv(f'/home/users/nus/e1124850/scratch/github/heritability_enrichment_ldsc/round5sqtl/tianchi_vat_{celltype}_gene_set_file.txt', index=False, header=False)





# baseldfiles = glob('/home/users/nus/e1124850/e1124850/co_lab/tianchi/1000G_Phase3_EAS_baselineLD_v2.2_ldscores/baselineLD.*.l2.ldscore.gz')

# total_SNP_df = pd.DataFrame()
# for ldfile in baseldfiles:
#     ld_df = pd.read_csv(ldfile, sep='\t', usecols=['SNP'])
#     total_SNP_df = pd.concat([total_SNP_df, ld_df['SNP']], axis=0)


# total_SNP_df.to_csv('/home/users/nus/e1124850/e1124850/co_lab/tongyihan/baselineLD_v2.2.snps.txt',header=False, index=False)

gene_coord_file_path = '/home/users/nus/e1124850/scratch/github/heritability_enrichment_ldsc/round5sqtl/gene_coord_file_hg37.txt'

##################
# sat_celltypes = sat_sqtl_df['celltype'].unique()
sat_celltypes = ['Adipocyte', 'Committed_preadipocyte', 'CEC', 'VEC', 'AEC', 'PVM',
       'SMC', 'Early_preadipocyte', 'LAM', 'Areg', 'Pericyte']

for celltype in sat_celltypes:
    sat_sqtl_sig_gene_path = f'/home/users/nus/e1124850/scratch/github/heritability_enrichment_ldsc/round5sqtl/tianchi_sat_{celltype}_gene_set_file.txt'
    for chrid in range(1, 23):
        cmd = (
            f'python /home/users/nus/e1124850/scratch/github/ldsc/make_annot.py '
            f'--gene-set-file {sat_sqtl_sig_gene_path} '
            f'--gene-coord-file {gene_coord_file_path} '
            f'--windowsize 100000 '
        f'--bimfile /home/users/nus/e1124850/e1124850/co_lab/tongyihan/round1/1000G_Phase3_EAS_plinkfiles/1000G.EAS.QC.{chrid}.bim '
        f'--annot-file /home/users/nus/e1124850/e1124850/co_lab/tianchi/gene_annot_file_grch37/tianchi_sat_{celltype}_sqtlsig_gene_annot_file.{chrid}.annot.gz'
    )
        os.system(cmd)


# vat_celltypes = vat_sqtl_df['celltype'].unique()

vat_celltypes = ['Mesothelial', 'Adipocyte', 'Committed_preadipocyte', 'PVM', 'CEC',
       'AEC', 'SMC', 'VEC', 'Early_preadipocyte', 'TIM4+_ATM', 'LEC',
       'Pericyte', 'IGFBP2+_cell', 'CD8+_T']

for celltype in vat_celltypes:
    vat_sqtl_sig_gene_path = f'/home/users/nus/e1124850/scratch/github/heritability_enrichment_ldsc/round5sqtl/tianchi_vat_{celltype}_gene_set_file.txt'
    for chrid in range(1, 23):
        cmd = (
            f'python /home/users/nus/e1124850/scratch/github/ldsc/make_annot.py '
            f'--gene-set-file {vat_sqtl_sig_gene_path} '
            f'--gene-coord-file {gene_coord_file_path} '
            f'--windowsize 100000 '
            f'--bimfile /home/users/nus/e1124850/e1124850/co_lab/tongyihan/round1/1000G_Phase3_EAS_plinkfiles/1000G.EAS.QC.{chrid}.bim '
            f'--annot-file /home/users/nus/e1124850/e1124850/co_lab/tianchi/gene_annot_file_grch37/tianchi_vat_{celltype}_sqtlsig_gene_annot_file.{chrid}.annot.gz'
        )
        os.system(cmd)
