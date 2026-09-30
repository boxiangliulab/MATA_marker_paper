import pathlib
import scanpy as sc
import pandas as pd
import os
from scipy import sparse
from scipy.io import mmwrite

# Project root is 8 levels up from this script
PROJECT_ROOT = pathlib.Path(__file__).resolve().parents[8]

adata_fp = str(PROJECT_ROOT / "data/08.integration_snRNA_object/01.freeze_version/SATA_integrated_obj_83_lib_freeze.h5ad")
adata = sc.read(adata_fp)
scvi_paths = {
    "SAT": "/home/project/11003054/e1353486/project/adipose/slingshot/02.integration/20260203225032/SATA_integrated_obj_83_lib_scvi_embed.csv",
    "VAT": "/home/project/11003054/e1353486/project/adipose/slingshot/02.integration/20260203154141/SATA_integrated_obj_83_lib_scvi_embed.csv"
}

out_base = str(PROJECT_ROOT / "data/09.eQTL/03.Saige_QTL/02.output/11.cell_state_eQTL/02.adipogenesis")

# ensure barcode index is string
adata.obs_names = adata.obs_names.astype(str)

# =====================
# Loop over depot
# =====================

for depot, scvi_file in scvi_paths.items():
    print(f"\nProcessing {depot}")
    # output directory
    outdir = os.path.join(out_base, depot)
    os.makedirs(outdir, exist_ok=True)
    # load scVI embedding
    scvi = pd.read_csv(scvi_file)
    barcodes = scvi.iloc[:, 0].astype(str).values
    print("Cells in scVI:", len(barcodes))
    # subset adata
    adata_sub = adata[barcodes].copy()
    print("Cells in adata subset:", adata_sub.n_obs)
    # counts matrix
    counts = adata_sub.layers["counts"]
    # ensure sparse
    if not sparse.issparse(counts):
        counts = sparse.csr_matrix(counts)
    # Save matrix
    # =====================
    mmwrite(os.path.join(outdir, "counts.mtx"), counts.T)
    # genes
    genes = pd.DataFrame({
        "gene_id": adata_sub.var_names,
        "gene_name": adata_sub.var_names
    })
    genes.to_csv(os.path.join(outdir, "genes.tsv"),
                 sep="\t",
                 index=False,
                 header=False)
    # barcodes
    pd.Series(adata_sub.obs_names).to_csv(
        os.path.join(outdir, "barcodes.tsv"),
        index=False,
        header=False
    )
    print(f"{depot} done ✓")

print("\nAll finished 🚀")
