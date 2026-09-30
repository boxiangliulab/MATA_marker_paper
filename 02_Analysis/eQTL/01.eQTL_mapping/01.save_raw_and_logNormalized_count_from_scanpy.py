# Date: 2025.4.1
# Author: Yihan Tong
# Usage: To save raw count and logNomalized data fron scanpy object


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
from matplotlib import gridspec
import numpy as np


####################### Save the count matrix of each cell type #############
# Project root is 8 levels up from this script
PROJECT_ROOT = pathlib.Path(__file__).resolve().parents[8]

adata_filtered = sc.read(PROJECT_ROOT / "data/08.integration_snRNA_object/01.freeze_version/SATA_integrated_obj_83_lib_freeze.h5ad")
marker_gene = pd.read_csv(PROJECT_ROOT / "data/08.integration_snRNA_object/00.annotation/02.annotation_level2/SATA_Annotation_level2_v3.csv")
cell_types = adata_filtered.obs['Annotation_level2'].dropna().unique()
outdir = PROJECT_ROOT / "data/08.integration_snRNA_object/00.annotation/02.annotation_level2/00.coordinate_per_cell_type/04.freeze_version"


## Save raw count and X by depot and cell type
def save_count_per_depot_celltype(adata_filtered, depot):
    """
    Save raw count data and normalized data for each depot and cell type
    Parameters:
    - adata_filtered: adata contains raw count layer and normalized and Logarithmize data
    - Depot: SAT or VAT
    Output:
    - raw count data and normalized data for each depot and cell type
    """
    output_dir = "/home/users/nus/e1101943/scratch/project/00.adipose/data/07.full_scale_snRNA/08.integration/all/scvi/03.count_matrix/"
    subset = adata_filtered[adata_filtered.obs['Depot'] == depot, :]
    cell_types = subset.obs['Annotation_level2'].dropna().unique()
    for cell_type in cell_types:
        adata_subset = subset[subset.obs['Annotation_level2'] == cell_type].copy()
        # raw count
        expr_df = pd.DataFrame(
        data=adata_subset.layers['counts'].todense(),
        index=adata_subset.obs.index,
        columns=adata_subset.var_names)
        output_path = os.path.join(output_dir, depot, "01.raw_count", f"{cell_type}.raw_count.tsv")
        expr_df.to_csv(output_path,sep='\t')
        # lognormalized data
        # expr_df = pd.DataFrame(
        # data=adata_subset.X.todense(),
        # index=adata_subset.obs.index,
        # columns=adata_filtered.var_names)
        # output_path = os.path.join(output_dir, depot, "02.logNormalize_data", f"{cell_type}.logNormalize_count.tsv")
        # expr_df.to_csv(output_path,sep='\t')
        print(f"Finished processing {cell_type}\n")

# save SAT
save_count_per_depot_celltype(adata_filtered,depot = "SAT")

# save VAT
save_count_per_depot_celltype(adata_filtered,depot = "VAT")
