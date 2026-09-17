#!/usr/bin/env python3
"""KNN meta-cell assignment within donor and depot.

Single nuclei are too sparse to quantify splicing, so nuclei of the same
cell type are pooled into "meta-cells". Meta-cells are built separately for
every donor and depot, so that each meta-cell corresponds to exactly one
donor / depot / cell-type combination and can be treated as one sample in
the sQTL model.

1. Take the scVI latent space (obsm["X_scVI"]) of the nuclei of one donor.
2. Greedily seed a meta-cell at each unassigned nucleus and fill it with its
   nearest unassigned neighbours until it holds `--meta_size` nuclei.
   Groups that cannot reach `--min_final_size` nuclei are discarded.
3. Nuclei left unassigned at the end are attached to the meta-cell whose
   centroid is closest in the latent space, so that no nucleus is lost.
"""

import argparse
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from sklearn.metrics import pairwise_distances
from sklearn.neighbors import NearestNeighbors

# --- Parse command-line arguments -------------------------------------------
parser = argparse.ArgumentParser(
    description="Construct meta-cells using KNN + centroid fallback")
parser.add_argument("--input_h5ad", type=str, required=True,
                    help="Path to the input h5ad file")
parser.add_argument("--cell_type", type=str, required=True,
                    help="Cell type name (e.g. Adipocyte, SMC)")
parser.add_argument("--depot", type=str, required=True,
                    help="Depot to process (SAT or VAT)")
parser.add_argument("--meta_size", type=int, default=10,
                    help="Target number of nuclei per meta-cell")
parser.add_argument("--min_final_size", type=int, default=6,
                    help="Minimum number of nuclei for a valid meta-cell")
parser.add_argument("--output_dir", type=str, required=True,
                    help="Directory to save output files")
args = parser.parse_args()

input_h5ad = args.input_h5ad
cell_type = args.cell_type
depot = args.depot
meta_size = args.meta_size
min_final_size = args.min_final_size

out_dir = args.output_dir
os.makedirs(out_dir, exist_ok=True)


def assign_leftover_by_centroid(X, meta_cell_labels, unassigned_indices):
    """Attach unassigned nuclei to the meta-cell with the nearest centroid.

    X                 : latent-space coordinates of all nuclei of this donor
    meta_cell_labels  : {meta-cell ID: [row indices of its nuclei]}
    unassigned_indices: row indices of nuclei not yet in any meta-cell
    """
    meta_ids = list(meta_cell_labels.keys())
    centroids = np.vstack([
        X[meta_cell_labels[mid]].mean(axis=0) for mid in meta_ids
    ])
    X_unassigned = X[unassigned_indices]
    dist = pairwise_distances(X_unassigned, centroids, metric="euclidean")
    for i, cell_idx in enumerate(unassigned_indices):
        nearest_meta = meta_ids[np.argmin(dist[i])]
        meta_cell_labels[nearest_meta].append(cell_idx)
    return meta_cell_labels


# --- Load data and subset to the requested depot -----------------------------
adata = sc.read_h5ad(input_h5ad)
adata.obs["meta_cell"] = ""
meta_info = []

if depot not in adata.obs["Depot"].unique():
    raise ValueError(f"{depot} not found in adata.obs['Depot']")

ad_depot = adata[adata.obs["Depot"] == depot].copy()

# --- Build meta-cells one donor at a time ------------------------------------
for donor in ad_depot.obs["donor_id"].unique():
    ad_donor = ad_depot[ad_depot.obs["donor_id"] == donor].copy()

    X = ad_donor.obsm["X_scVI"]
    cell_names = np.array(ad_donor.obs_names)

    # Donors with too few nuclei cannot form even one meta-cell.
    if len(cell_names) < min_final_size:
        continue

    # Nearest neighbours of every nucleus in the scVI latent space.
    nbrs = NearestNeighbors(
        n_neighbors=min(meta_size, len(cell_names)), algorithm="auto").fit(X)
    distances, indices = nbrs.kneighbors(X)

    assigned = set()
    group_list = []

    # Greedy seeding: walk over the nuclei and grow a meta-cell around each
    # nucleus that has not been assigned yet.
    for i in range(len(cell_names)):
        if cell_names[i] in assigned:
            continue
        neighbors = indices[i][1:]          # skip the nucleus itself
        group = [cell_names[i]]
        for j in neighbors:
            candidate = cell_names[j]
            if candidate not in assigned:
                group.append(candidate)
            if len(group) == meta_size:
                break
        if len(group) < min_final_size:     # too small to keep
            continue
        group_list.append(group)
        assigned.update(group)

    unassigned = list(set(cell_names) - assigned)

    if len(unassigned) > 0:
        # Leftover nuclei join the closest existing meta-cell.
        meta_cell_labels = {
            f"{donor}_{depot}_{cell_type}_MC{idx + 1:02d}":
                [np.where(cell_names == cell)[0][0] for cell in group]
            for idx, group in enumerate(group_list)
        }
        unassigned_indices = [np.where(cell_names == cell)[0][0]
                              for cell in unassigned]
        if len(meta_cell_labels) > 0:
            meta_cell_labels = assign_leftover_by_centroid(
                X, meta_cell_labels, unassigned_indices)
            for meta_id, idx_list in meta_cell_labels.items():
                cluster_id = int(meta_id.split("_MC")[-1])
                for idx in idx_list:
                    cell = cell_names[idx]
                    adata.obs.at[cell, "meta_cell"] = meta_id
                    meta_info.append({
                        "cell": cell,
                        "donor": donor,
                        "depot": depot,
                        "meta_cell": meta_id,
                        "meta_cluster": cluster_id,
                    })
    else:
        # Every nucleus was assigned during the greedy pass.
        for idx, group in enumerate(group_list):
            meta_id = f"{donor}_{depot}_{cell_type}_MC{idx + 1:02d}"
            for cell in group:
                adata.obs.at[cell, "meta_cell"] = meta_id
                meta_info.append({
                    "cell": cell,
                    "donor": donor,
                    "depot": depot,
                    "meta_cell": meta_id,
                    "meta_cluster": idx + 1,
                })

    print(f"{donor}-{depot}: processed {len(cell_names)} cells")

# --- Save the nucleus -> meta-cell map ---------------------------------------
meta_df = pd.DataFrame(meta_info)
meta_df.to_csv(
    f"{out_dir}/{cell_type}_{depot}_meta_cell_knn_fixed_{meta_size}_centroid.csv",
    index=False)

# --- QC plot 1: distribution of meta-cell sizes ------------------------------
meta_counts = meta_df.groupby(["depot", "meta_cell"]).size().reset_index(
    name="n_cells")

min_count = meta_counts["n_cells"].min()
max_count = meta_counts["n_cells"].max()

plt.figure(figsize=(8, 4))
sns.histplot(
    data=meta_counts,
    x="n_cells",
    element="step",
    bins=range(min_count, max_count + 2),
    common_norm=False,
)
plt.title(f"{cell_type} ({depot}): Meta-cell Size Distribution")
plt.xlabel("Number of cells per meta-cell")
plt.ylabel("Frequency")
plt.tight_layout()
plt.savefig(
    f"{out_dir}/{cell_type}_{depot}_meta_cell_size_distribution_{meta_size}.png")
plt.savefig(
    f"{out_dir}/{cell_type}_{depot}_meta_cell_size_distribution_{meta_size}.pdf")
plt.close()

# --- QC plot 2: number of meta-cells contributed by each donor ---------------
meta_donor_counts = meta_df.groupby(["depot", "donor"]).meta_cell.nunique(
).reset_index(name="n_meta_cells")

plt.figure(figsize=(8, 4))
sns.histplot(
    data=meta_donor_counts,
    x="n_meta_cells",
    element="step",
    bins=np.logspace(np.log10(meta_donor_counts["n_meta_cells"].min()),
                     np.log10(meta_donor_counts["n_meta_cells"].max() + 1),
                     num=20),
    common_norm=False,
)
plt.xscale("log")
plt.title(f"{cell_type} ({depot}): Number of Meta-cells per Donor")
plt.xlabel("Number of meta-cells per donor")
plt.ylabel("Frequency")
plt.tight_layout()
plt.savefig(
    f"{out_dir}/{cell_type}_{depot}_meta_cell_count_per_donor_{meta_size}.png")
plt.savefig(
    f"{out_dir}/{cell_type}_{depot}_meta_cell_count_per_donor_{meta_size}.pdf")
plt.close()

print(f"Meta-cell assignment written to {out_dir}")
