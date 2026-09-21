import pandas as pd
import os

# 每个 cell type 取 ES 最高的前 10% 基因,其余基因的 ES 置 0,全部 cell type 合成一个矩阵。
#
# 为什么是"置 0"而不是"删掉基因行":
#   CELLECT 生成注释时,矩阵里不存在的基因一律按 ES=0 处理
#   (make_annot_from_geneset_all_chr_snake.py 的 get_max_ES_vals),
#   所以两种写法产生的注释完全等价。
#
# 为什么必须合成一个矩阵:
#   LDSC 在 n_annot==1 时会把 ldscore 的列名直接写成 "L2",把注释名丢掉
#   (ldsc.py: if n_annot == 1: ldscore_colnames = [col_prefix+scale_suffix]),
#   于是 CELLECT 的 split_ldscores_snake.py 做 re.sub(r"L2$","",x) 得到空串,
#   报 KeyError: "None of [Index([''])] are in the [columns]"。
#   一个矩阵放多列注释即可绕开,同时 LD score 只需算一次而不是每个 cell type 算一遍,
#   all-genes 背景注释也统一为全部基因,cell type 之间可比。

filtered_top_genes_percent = 0.1
outdir = '/home/users/nus/e1124850/e1124850/co_lab/tongyihan/round4'

cellex_files = {
    'SAT': '/data/projects/11003054/e1101943/project/01.adipose_sqtl/data/33.CELLEX/cellex_Annotation_level2.SAT_filtered.ensembl.esmu.csv.gz',
    'VAT': '/data/projects/11003054/e1101943/project/01.adipose_sqtl/data/33.CELLEX/cellex_Annotation_level2.VAT_filtered.ensembl.esmu.csv.gz',
}

for tissue, cellex_file in cellex_files.items():
    df = pd.read_csv(cellex_file)
    # CELLECT 的 id/annotation 名只允许字母数字和单个下划线/连字符
    df = df.rename(columns={c: c.replace('+', 'pos') for c in df.columns})
    df = df.set_index('gene')

    n_top = int(len(df) * filtered_top_genes_percent)
    top_df = pd.DataFrame(0.0, index=df.index, columns=df.columns)
    for celltype in df.columns:
        top_genes = df[celltype].sort_values(ascending=False).index[:n_top]
        top_df.loc[top_genes, celltype] = df.loc[top_genes, celltype]

    output = os.path.join(outdir, f'{tissue}_top_{int(filtered_top_genes_percent*100)}percent.csv.gz')
    top_df.to_csv(output, compression='gzip')
    print(output, top_df.shape, 'top genes per cell type:', n_top)
