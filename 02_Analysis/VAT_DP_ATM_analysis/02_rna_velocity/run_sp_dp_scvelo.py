#!/usr/bin/env python
# PURPOSE: Fit the scVelo dynamical model for SP and DP ATMs and infer velocity
#          pseudotime and direction on the VAT macrophage embedding.
# INPUT: Annotated macrophage H5AD plus per-library velocyto loom files.
# OUTPUT: Processed H5AD, fate probabilities, and CellRank state summaries.
# USAGE: Set the variables in config/paths.example.env and run this script in
#        the Linux/HPC environment used to generate the velocyto loom files.
"""
Infer directionality between TIM4+ (SP) and TIM4+CD11c+ (DP) ATMs from
spliced and unspliced counts. CellRank complements scVelo by estimating
data-driven initial states, terminal states, and fate probabilities.

Loom CellIDs use ``sample_lane:barcodex`` whereas H5AD observation names use
``barcode_sample_lane``. The conversion below was verified against the source
data before the original analysis.
"""
import os
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"  # Avoid read-only NFS locking errors.
import glob
import numpy as np
# Compatibility shim for CellRank 2.0.7 with NumPy 2.x.
import numpy.testing as _npt
_orig_aae = _npt.assert_array_equal
def _aae_compat(*args, **kw):
    if "x" in kw: kw["actual"] = kw.pop("x")
    if "y" in kw: kw["desired"] = kw.pop("y")
    return _orig_aae(*args, **kw)
_npt.assert_array_equal = _aae_compat
np.testing.assert_array_equal = _aae_compat
# --- end shim ---
import pandas as pd
import scanpy as sc
import scvelo as scv
import anndata as ad

scv.settings.verbosity = 3
scv.settings.presenter_view = True

VELO_DIR = os.path.join(os.environ["SECTION4_OUTPUT"], "rna_velocity")
LOOM_DIR = os.environ["SECTION4_VELOCYTO_LOOMS"]
RES_DIR  = os.path.join(VELO_DIR, "results")
os.makedirs(RES_DIR, exist_ok=True)

H5AD_ANNOT = os.environ["SECTION4_MACROPHAGE_H5AD"]

ANNOT_COL = "Annotation_level2"
LIB_COL   = "Library"                 # Example: B10_L1; matches loom filenames.
SP_LABEL  = "TIM4+_ATM"               # single positive
DP_LABEL  = "TIM4+_CD11c+_ATM"        # double positive

# Reuse the merged checkpoint to avoid re-reading all loom files.
CKPT = os.path.join(RES_DIR, "sp_dp_merged_raw.h5ad")

if os.path.exists(CKPT):
    print(f"[checkpoint] Loading {CKPT}; skipping loom import")
    adata = sc.read_h5ad(CKPT)
    print("    loaded:", adata.shape, "| layers:", list(adata.layers.keys()))
else:
    # 1. Load the annotated macrophage object and retain SP and DP ATMs.
    print("[1] Loading annotated H5AD")
    adata = sc.read_h5ad(H5AD_ANNOT)
    print("    all macrophages:", adata.shape)
    adata = adata[adata.obs[ANNOT_COL].isin([SP_LABEL, DP_LABEL])].copy()
    adata.obs[ANNOT_COL] = adata.obs[ANNOT_COL].astype("category")
    print("    SP/DP subset:", adata.shape)
    print(adata.obs[ANNOT_COL].value_counts().to_string())

    needed = set(adata.obs_names)
    libs_needed = set(adata.obs[LIB_COL].astype(str).unique())
    print(f"    required cells: {len(needed)}; libraries: {len(libs_needed)}")

    # 2. Read required cells from each loom file, then concatenate them.
    def loom_to_keys(cellids):
        # "B10_L1:AAAC...x" -> "AAAC..._B10_L1"
        out = []
        for c in cellids:
            samp, bc = c.split(":")
            out.append(bc[:-1] + "_" + samp if bc.endswith("x") else bc + "_" + samp)
        return np.array(out)

    loom_files = sorted(glob.glob(os.path.join(LOOM_DIR, "*.loom")))
    print(f"[2] Scanning {len(loom_files)} loom files")
    pieces = []
    for f in loom_files:
        samp = os.path.basename(f).replace(".loom", "")
        if samp not in libs_needed:
            continue
        ld = sc.read_loom(f)
        ld.var_names_make_unique()
        ld.obs_names = loom_to_keys(ld.obs_names.values)
        ld.obs_names_make_unique()
        keep = ld.obs_names.isin(needed)
        n = int(keep.sum())
        if n == 0:
            continue
        pieces.append(ld[keep].copy())
        print(f"    {samp}: retained {n} cells")
        del ld

    ldata = pieces[0] if len(pieces) == 1 else ad.concat(pieces, join="outer", merge="same")
    ldata.obs_names_make_unique()
    print("    combined loom:", ldata.shape, "| layers:", list(ldata.layers.keys()))

    # 3. Add spliced and unspliced layers to the annotated object.
    print("[3] merge spliced/unspliced into annotated adata")
    adata = scv.utils.merge(adata, ldata)
    print("    merged:", adata.shape, "| layers:", list(adata.layers.keys()))
    print(adata.obs[ANNOT_COL].value_counts().to_string())
    adata.write(CKPT)
    print(f"    saved checkpoint: {CKPT}")

# ============================================================
# 4. Preprocess counts and compute moments. Scanpy functions are used explicitly
# because the original scVelo build had an n_top_genes forwarding issue.
# ============================================================
print("[4] filter genes + normalize + log1p + HVG + moments")
scv.pp.filter_genes(adata, min_shared_counts=20)
scv.pp.normalize_per_cell(adata)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, subset=True)
print("    HVG subset:", adata.shape)
sc.pp.pca(adata, n_comps=50)
scv.pp.moments(adata, n_pcs=50, n_neighbors=30)

# ============================================================
# 5. Fit the dynamical velocity model and estimate latent time.
# ============================================================
print("[5] recover_dynamics + velocity(dynamical)")
scv.tl.recover_dynamics(adata, n_jobs=8)
scv.tl.velocity(adata, mode="dynamical")
scv.tl.velocity_graph(adata, n_jobs=8)
scv.tl.latent_time(adata)

# ============================================================
# 6. Compute velocity pseudotime. Figure rendering is intentionally excluded
# from this analysis-only release.
# ============================================================
print("[6] Computing velocity pseudotime")
scv.tl.velocity_pseudotime(adata)

# ============================================================
# 7. Estimate data-driven macrostates and fate probabilities with CellRank.
# ============================================================
print("[7] CellRank fate probabilities")
import cellrank as cr
from cellrank.kernels import VelocityKernel, ConnectivityKernel

vk = VelocityKernel(adata).compute_transition_matrix()
ck = ConnectivityKernel(adata).compute_transition_matrix()
combined = 0.8 * vk + 0.2 * ck

g = cr.estimators.GPCCA(combined)
g.compute_macrostates(n_states=2, cluster_key=ANNOT_COL)

# A higher coarse-grained stationary probability supports a terminal state;
# initial and terminal states are inferred without manually assigning direction.
coarse_stat = getattr(g, "coarse_stationary_distribution", None)
coarse_T    = getattr(g, "coarse_T", None)
print("[7a] Coarse-grained stationary distribution:")
print(coarse_stat)
print("[7b] coarse-grained transition matrix:")
print(coarse_T)

# Infer initial and terminal states directly from the transition structure.
g.predict_terminal_states()
g.predict_initial_states(allow_overlap=True)
g.compute_fate_probabilities()

with open(os.path.join(RES_DIR, "cellrank_states.txt"), "w") as fh:
    fh.write("=== macrostates ===\n" + str(g.macrostates.value_counts()) + "\n\n")
    fh.write("=== coarse stationary distribution (higher supports terminal state) ===\n")
    fh.write(str(coarse_stat) + "\n\n")
    fh.write("=== coarse transition matrix ===\n" + str(coarse_T) + "\n\n")
    fh.write("=== terminal states (data-driven) ===\n" + str(g.terminal_states.cat.categories.tolist()) + "\n")
    fh.write("=== initial states  (data-driven) ===\n"  + str(g.initial_states.cat.categories.tolist()) + "\n")

fp = g.fate_probabilities
fate_df = pd.DataFrame(fp, index=adata.obs_names, columns=[str(c) for c in fp.names])
fate_df[ANNOT_COL]  = adata.obs[ANNOT_COL].values
fate_df["latent_time"] = adata.obs["latent_time"].values
fate_df["velocity_pseudotime"] = adata.obs["velocity_pseudotime"].values
fate_df.to_csv(os.path.join(RES_DIR, "sp_dp_fate_probabilities.csv"))

print("Mean time by label (lower values support an earlier state):")
print(fate_df.groupby(ANNOT_COL)[["latent_time", "velocity_pseudotime"]].mean().to_string())

# ============================================================
# 8. Convert CellRank lineage objects to arrays before writing with AnnData 0.12.
# ============================================================
for k in list(adata.obsm.keys()):
    v = adata.obsm[k]
    if not isinstance(v, np.ndarray):
        try:
            adata.obsm[k] = np.asarray(v.X) if hasattr(v, "X") else np.asarray(v)
        except Exception:
            del adata.obsm[k]
try:
    adata.write(os.path.join(RES_DIR, "sp_dp_velocity.h5ad"))
    print("[8] Saved processed H5AD")
except Exception as e:
    print("[8] H5AD write failed; tabular outputs were already saved:", e)
