#!/usr/bin/env python3
"""Identify marker genes across macrophage subtypes.

Input: annotated H5AD with Depot, Annotation_level1, and Annotation_level2.
Output: complete Scanpy one-versus-rest marker statistics as CSV or TSV.
By default, SAT and VAT macrophages are analysed together, matching the DEG
input used for the macrophage-subtype GSEA. Use ``--depot VAT`` for a
VAT-specific marker analysis.
No diagnostic or manuscript figures are produced.

Usage:
  python macrophage_marker_deg.py --h5ad INPUT.h5ad --output DEG_lv2_subtypes.csv
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import scanpy as sc


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5ad", required=True, help="Annotated single-cell H5AD")
    parser.add_argument("--output", required=True, help="Output marker TSV")
    parser.add_argument("--depot", default="all", help="Depot to retain, or 'all'")
    parser.add_argument("--major-label", default="Macrophage")
    parser.add_argument("--method", default="wilcoxon", choices=["wilcoxon", "t-test"])
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    adata = sc.read_h5ad(args.h5ad)
    required = {"Depot", "Annotation_level1", "Annotation_level2"}
    missing = required.difference(adata.obs.columns)
    if missing:
        raise ValueError(f"Missing observation columns: {sorted(missing)}")

    keep = adata.obs["Annotation_level1"] == args.major_label
    if args.depot.lower() != "all":
        keep &= adata.obs["Depot"] == args.depot
    adata = adata[keep].copy()
    if adata.n_obs == 0:
        raise ValueError("No cells remain after filtering")

    adata.obs["Annotation_level2"] = (
        adata.obs["Annotation_level2"].astype("category").cat.remove_unused_categories()
    )
    sc.tl.rank_genes_groups(
        adata,
        groupby="Annotation_level2",
        method=args.method,
        pts=True,
    )

    tables = []
    for subtype in adata.obs["Annotation_level2"].cat.categories:
        table = sc.get.rank_genes_groups_df(adata, group=subtype)
        table.insert(0, "group", subtype)
        tables.append(table)

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    result = pd.concat(tables, ignore_index=True)
    if output.suffix.lower() == ".csv":
        result.to_csv(output, index=False)
    else:
        result.to_csv(output, sep="\t", index=False)
    print(f"Saved {output}")


if __name__ == "__main__":
    main()
