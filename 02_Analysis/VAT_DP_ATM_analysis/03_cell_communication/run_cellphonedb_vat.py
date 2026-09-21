#!/usr/bin/env python3
"""Run CellPhoneDB v5 on annotated VAT cells through OmicVerse.

Inputs are provided through SECTION4_VAT_H5AD and
SECTION4_CELLPHONEDB_REFERENCE. Results are written below SECTION4_OUTPUT.
The script performs the statistical interaction analysis only; visualization
and figure assembly are intentionally excluded from this release.

Usage: set the variables in config/paths.example.env and run this script.
"""

import os
import multiprocessing as mp

def main():
    # Fork avoids recursive worker initialization in the Linux HPC environment.
    try:
        mp.set_start_method("fork", force=True)
    except RuntimeError:
        # start method already set
        pass

    import scanpy as sc
    import omicverse as ov


    # Output paths
    TAG = os.environ.get("SECTION4_CELLPHONEDB_TAG", "vat")
    base_out = os.path.join(os.environ["SECTION4_OUTPUT"], f"{TAG}_cellphonedb")

    cpdb_raw_dir = os.path.join(base_out, "cpdb_raw")
    for d in [base_out, cpdb_raw_dir]:
        os.makedirs(d, exist_ok=True)

    pkl_path  = os.path.join(base_out, f"cpdb_vat_{TAG}.pkl")
    h5ad_path = os.path.join(base_out, f"cpdb_vat_{TAG}.h5ad")

    # Load annotated VAT cells.
    adata_path = os.environ["SECTION4_VAT_H5AD"]
    adata = ov.read(adata_path)

    # CellPhoneDB permutation analysis.
    cpdb_results, adata_cpdb = ov.single.run_cellphonedb_v5(
        adata,
        cpdb_file_path=os.environ["SECTION4_CELLPHONEDB_REFERENCE"],
        celltype_key="Annotation_level2",
        min_cell_fraction=0.01,
        min_genes=200,
        min_cells=3,
        iterations=1000,
        threshold=0.1,
        pvalue=0.05,
        threads=5,
        output_dir=cpdb_raw_dir,
        cleanup_temp=True
    )

    # Save raw and serialized results for downstream extraction.
    ov.utils.save(cpdb_results, pkl_path)
    adata_cpdb.write(h5ad_path, compression="gzip")

    print(f"Saved CellPhoneDB results to {base_out}")

if __name__ == "__main__":
    main()
