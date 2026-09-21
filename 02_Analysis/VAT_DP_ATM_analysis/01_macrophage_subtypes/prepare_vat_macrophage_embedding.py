#!/usr/bin/env python3
"""Recompute the VAT macrophage neighbor graph and UMAP coordinates.

The script writes an analysis H5AD and a cell-level coordinate/expression TSV.
It does not generate a figure.

Usage:
  python prepare_vat_macrophage_embedding.py --h5ad INPUT.h5ad --output OUTPUT_DIR
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5ad", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--embedding", default="X_scVI")
    parser.add_argument("--markers", default="TIMD4,MARCO,ITGAX,IL1B")
    args = parser.parse_args()

    output = Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    adata = sc.read_h5ad(args.h5ad)
    required = {"Depot", "Annotation_level1", "Annotation_level2"}
    missing = required.difference(adata.obs.columns)
    if missing:
        raise ValueError(f"Missing observation columns: {sorted(missing)}")
    if args.embedding not in adata.obsm:
        raise ValueError(f"Embedding {args.embedding!r} was not found in .obsm")

    keep = (adata.obs["Depot"].astype(str) == "VAT") & (
        adata.obs["Annotation_level1"].astype(str) == "Macrophage"
    )
    adata = adata[keep].copy()
    sc.pp.neighbors(adata, use_rep=args.embedding)
    sc.tl.umap(adata, random_state=0)
    adata.write(output / "vat_macrophage_embedding.h5ad", compression="gzip")

    markers = [x for x in args.markers.split(",") if x in adata.var_names]
    table = adata.obs.copy()
    table["UMAP_1"] = adata.obsm["X_umap"][:, 0]
    table["UMAP_2"] = adata.obsm["X_umap"][:, 1]
    if markers:
        expression = adata[:, markers].X
        if sp.issparse(expression):
            expression = expression.toarray()
        expression = np.asarray(expression, dtype=np.float32)
        for index, gene in enumerate(markers):
            table[f"{gene}_log1p"] = expression[:, index]
            sd = expression[:, index].std(ddof=0)
            table[f"{gene}_scaled"] = (
                0 if sd == 0 else (expression[:, index] - expression[:, index].mean()) / sd
            )
    table.to_csv(output / "vat_macrophage_embedding.tsv", sep="\t", index_label="cell_id")


if __name__ == "__main__":
    main()
