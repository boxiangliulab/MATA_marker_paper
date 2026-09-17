#!/usr/bin/env python3
"""Number of nuclei per donor and depot for one cell type.

Summarise how many nuclei each donor contributes to a given cell type in
each depot (SAT / VAT) and derive the target meta-cell size used in step 02.

The target size is ceil(median nuclei per donor / 4), i.e. the median donor
is split into roughly four meta-cells. Donors with zero nuclei are excluded
from the median so that absent donors do not shrink the meta-cell size.
"""

import argparse
import os

import numpy as np
import pandas as pd
import scanpy as sc

# --- Parse arguments ---------------------------------------------------------
parser = argparse.ArgumentParser(
    description="Summarise nuclei per donor and depot for one cell type")
parser.add_argument("--input_h5ad", type=str, required=True,
                    help="Path to the input h5ad file")
parser.add_argument("--cell_type", type=str, required=True,
                    help="Cell type name (e.g. Adipocyte, SMC)")
parser.add_argument("--output_dir", type=str, required=True,
                    help="Directory to save the output tables")
args = parser.parse_args()

out_dir = args.output_dir
os.makedirs(out_dir, exist_ok=True)

# --- Load data ---------------------------------------------------------------
adata = sc.read_h5ad(args.input_h5ad)

# --- Nuclei per donor and depot ---------------------------------------------
cell_counts = (
    adata.obs.groupby(["donor_id", "Depot"]).size().reset_index(name="n_cells")
)

# Order donors by total contribution (largest first) for readability.
donor_totals = cell_counts.groupby("donor_id")["n_cells"].sum().reset_index()
donor_order = donor_totals.sort_values("n_cells", ascending=False)["donor_id"]
cell_counts["donor_id"] = pd.Categorical(
    cell_counts["donor_id"], categories=donor_order, ordered=True)

cell_count_table_path = os.path.join(
    out_dir, f"{args.cell_type}_donor_depot_cell_counts.tsv")
cell_counts.to_csv(cell_count_table_path, sep="\t", index=False)

# --- Median per depot and the derived meta-cell size -------------------------
# Donors that contribute no nuclei to a depot are excluded from the median.
median_table = (
    cell_counts[cell_counts["n_cells"] > 0]
    .groupby("Depot")["n_cells"]
    .median()
    .reset_index()
    .rename(columns={"n_cells": "median_n_cells_per_donor"})
)
# Target meta-cell size: roughly four meta-cells for the median donor.
median_table["median_div4_ceil"] = np.ceil(
    median_table["median_n_cells_per_donor"] / 4).astype(int)

median_table_path = os.path.join(
    out_dir, f"{args.cell_type}_median_cells_per_donor.tsv")
median_table.to_csv(median_table_path, sep="\t", index=False)

print(f"Cell counts written to : {cell_count_table_path}")
print(f"Median table written to: {median_table_path}")
