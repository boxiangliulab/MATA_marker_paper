import pathlib
import scanpy as sc
import scipy.io
import scipy.sparse as sp
import pandas as pd
import os
import gzip
import shutil

# Project root is 7 levels up from this script
PROJECT_ROOT = pathlib.Path(__file__).resolve().parents[7]

H5AD      = str(PROJECT_ROOT / "data/08.integration_snRNA_object/01.freeze_version/SATA_integrated_obj_83_lib_freeze.h5ad")
OUT_BASE  = str(PROJECT_ROOT / "data/35.vascular/SMC_counts")
ASIAN_RACES = ["Chinese", "Malay", "Indian"]
CELL_TYPE   = "SMC"

# ---------- load --------------------------------------------------------------
print("Loading h5ad ...")
adata = sc.read_h5ad(H5AD)
print(f"  Total cells: {adata.n_obs}, genes: {adata.n_vars}")

# ---------- filter ------------------------------------------------------------
mask = (
    adata.obs["Race"].isin(ASIAN_RACES) &
    (adata.obs["Annotation_level2"] == CELL_TYPE)
)
adata_sub = adata[mask].copy()
print(f"  After filter (Asian + {CELL_TYPE}): {adata_sub.n_obs} cells")
print(adata_sub.obs["Depot"].value_counts())

# ---------- raw counts --------------------------------------------------------
if "counts" in adata_sub.layers:
    mat   = adata_sub.layers["counts"]
    genes = adata_sub.var_names
    print("  Using layers['counts']")
elif adata_sub.raw is not None:
    raw   = adata_sub.raw.to_adata()
    raw   = raw[adata_sub.obs_names].copy()
    mat   = raw.X
    genes = raw.var_names
    print("  Using adata.raw")
else:
    mat   = adata_sub.X
    genes = adata_sub.var_names
    print("  Using adata.X (no raw/counts layer found)")

# ---------- export per depot --------------------------------------------------
def save_10x(cells_mask, depot_name, count_mat, cell_obs, gene_names):
    out_dir = os.path.join(OUT_BASE, depot_name)
    os.makedirs(out_dir, exist_ok=True)

    sub_mat = count_mat[cells_mask]
    sub_obs = cell_obs[cells_mask]

    if not sp.issparse(sub_mat):
        sub_mat = sp.csc_matrix(sub_mat)
    else:
        sub_mat = sub_mat.tocsc()

    # barcodes
    pd.Series(sub_obs.index).to_csv(
        os.path.join(out_dir, "barcodes.tsv.gz"),
        index=False, header=False, compression="gzip"
    )

    # features
    feat = pd.DataFrame({
        "gene_id":       gene_names,
        "gene_name":     gene_names,
        "feature_type":  "Gene Expression"
    })
    feat.to_csv(
        os.path.join(out_dir, "features.tsv.gz"),
        index=False, header=False, sep="\t", compression="gzip"
    )

    # matrix (genes x cells)
    mtx_path = os.path.join(out_dir, "matrix.mtx")
    scipy.io.mmwrite(mtx_path, sub_mat.T)
    with open(mtx_path, "rb") as f_in, \
         gzip.open(mtx_path + ".gz", "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
    os.remove(mtx_path)

    # metadata
    sub_obs.to_csv(os.path.join(out_dir, "metadata.tsv.gz"),
                   sep="\t", compression="gzip")

    print(f"  Saved {depot_name}: {sub_mat.shape[0]} cells -> {out_dir}")


for depot in ["SAT", "VAT"]:
    mask_depot = adata_sub.obs["Depot"] == depot
    save_10x(
        cells_mask = mask_depot.values,
        depot_name = depot,
        count_mat  = mat,
        cell_obs   = adata_sub.obs,
        gene_names = genes
    )

print("Done.")
