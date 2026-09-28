import os
import tempfile
import sys

import scanpy.external as sce
import scanpy as sc
import scvi
import seaborn as sns
import torch
import matplotlib.pyplot as plt
import pandas as pd
import argparse
import pathlib
import yaml
import jax
import time
import uuid
import anndata
from scipy.sparse import csr_matrix
from scib_metrics.benchmark import Benchmarker

###### Pre-settings
scvi.settings.seed = 0
print("Last run with scvi-tools version:", scvi.__version__)

sc.set_figure_params(figsize=(6, 6), frameon=False)
sns.set_theme()
torch.set_float32_matmul_precision("high")

# ###### Reading files
# matrix_path_list = pd.read_csv("./scvi/00.filtered_matrix_path.txt", header = None, delimiter = "\t")
# ### [1] Convert count matrix into a h5ad file
# #adatas = []
# #for i in range(len(matrix_path_list[0])):
# #    current_adata = sc.read(matrix_path_list[0][i],delimiter='\t',cache = False)
# #    current_adata = current_adata.T
# #    library_name = matrix_path_list[0][i].split("/")[8].split("_")[0] + "_" + matrix_path_list[0][i].split("/")[8].split("_")[1]
# #    print(f"Number of cells in {library_name}: {current_adata.n_obs}")
# #    current_adata.obs['Library'] = library_name
# #    print(current_adata)
# #    current_adata.write(f"./01.scvi/{library_name}.h5ad", compression="gzip")
# #    adatas.append(current_adata)

# #adata = anndata.concat(adatas)
# #adata.write("./scvi/MATA_integrated_obj_83_lib.h5ad", compression="gzip")
# library_names = []
# for i in range(len(matrix_path_list[0])):
#      library_name = matrix_path_list[0][i].split("/")[8].split("_")[0] + "_" + matrix_path_list[0][i].split("/")[8].split("_")[1]
#      library_names.append(library_name)
#      len(library_names)


# print("Finishing writing h5ads")

# adatas = []
# for j in range(len(library_names)):
#      current_adata = sc.read(os.path.join("./01.scvi", f"{library_names[j]}.h5ad"))
#      current_adata.X = csr_matrix(current_adata.X)
#      adatas.append(current_adata)

# adata = anndata.concat(adatas, join='outer')
# metadata_df = pd.read_csv("./metadata_after_libraryQC_b4_ctQC.tsv", sep="\t")
# metadata_df = metadata_df.set_index("barcode")
# metadata_df = metadata_df.loc[adata.obs_names.intersection(metadata_df.index)]
# adata.obs = adata.obs.join(metadata_df, how="left")
# #adata.write("./scvi/MATA_integrated_obj_83_lib.h5ad", compression="gzip")

adata = anndata.read_h5ad("./scvi/MATA_integrated_obj_83_lib.h5ad")

save_path = "./scvi/"
#batch_key = "Library"
batch_key = "donor_id"
now = time.strftime("%Y%m%d%H%M%S", time.localtime(time.time()))
save_path = os.path.join(save_path, f"{now}")
if not os.path.exists(save_path):
	os.makedirs(save_path)

sc.pp.filter_genes(adata, min_counts=3)
# adata = adata.raw.to_adata()
adata.layers["counts"] = adata.X.copy()  # preserve counts
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

sc.pp.highly_variable_genes(
	adata,
	n_top_genes = 3000,
	subset = True,
	layer = "counts",
	flavor = 'seurat_v3_paper',
	batch_key = batch_key)
########################## scVI ###########################
scvi.model.SCVI.setup_anndata(
		adata,
		layer="counts",
		batch_key = batch_key
		)

model = scvi.model.SCVI(adata, n_hidden = 128, n_latent = 30, n_layers = 2, gene_likelihood="nb")
model.train(max_epochs = 400, batch_size = 128)

model_dir = os.path.join(save_path, 'SCVI_model')
model.save(model_dir, overwrite=True)

model = scvi.model.SCVI.load(model_dir, adata=adata)

SCVI_LATENT_KEY = "X_scVI"

latent = model.get_latent_representation()
adata.obsm[SCVI_LATENT_KEY] = latent

sc.pp.neighbors(adata, use_rep=SCVI_LATENT_KEY)
sc.tl.umap(adata)

########### scvi: prepare output files #####
sc.pl.umap(adata,
			color=["Library"],
			frameon=False,
			size=20
			)
plt.savefig(os.path.join(save_path, f"Batch_cell_type_UMAP.pdf"), bbox_inches='tight', dpi = 300)

SCVI_CLUSTERS_KEY = f"leiden_scVI"

sc.tl.leiden(adata, key_added=SCVI_CLUSTERS_KEY, resolution=1)
sc.pl.umap(
            adata,
            color=[SCVI_CLUSTERS_KEY],
            frameon=False,
            size=20,
            legend_loc='on data',
        )
plt.savefig(os.path.join(save_path, f"Clustering_Leiden_UMAP.pdf"), bbox_inches='tight', dpi=300)

adata.write(os.path.join(save_path, "MATA_integrated_obj_83_lib.h5ad"), compression="gzip")

import pandas as pd
df = pd.DataFrame(data=adata.obsm["X_scVI"], index=adata.obs_names)
df.to_csv(os.path.join(save_path,"MATA_integrated_obj_83_lib_scvi_embed.csv"))

df_umap = pd.DataFrame(data=adata.obsm["X_umap"], index=adata.obs_names)
df_umap.to_csv(os.path.join(save_path,"MATA_integrated_obj_83_lib_scvi_UMAP.csv"))

plt.savefig(os.path.join(save_path, f"scAVNI_Clustering_Leiden_UMAP.pdf"), bbox_inches='tight', dpi=300)
cluster = pd.DataFrame(data=adata.obs["leiden_scVI"], index=adata.obs_names)
cluster.to_csv(os.path.join(save_path,"Clustering_Leiden_res1_UMAP.csv"))
print("Finishing saving scvi output")

