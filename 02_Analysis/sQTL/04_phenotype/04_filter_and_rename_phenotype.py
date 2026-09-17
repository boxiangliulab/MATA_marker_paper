#!/usr/bin/env python3
"""Restrict and rename the meta-cell phenotype matrix.

The meta-cell PSI matrix produced in step 03 is unfiltered and its columns
are named after the meta-cell files. This script

  1. keeps only the splice sites (rows) that passed the pseudobulk filters
     -- the reference file given with -r -- so that the same set of sites is
     tested in every sample, and
  2. renames the columns to "<donor>:<n>", i.e. the donor the meta-cell came
     from plus a running index. The sQTL model uses these names to link the
     repeated meta-cell measurements of a donor to that donor's genotype.

The original -> new column names are written to the mapping file (-m).
"""

import os
import re
from collections import defaultdict

import pandas as pd


def read_and_fix_header(input_file):
    """Read a phenotype matrix whose header has one field fewer than the rows.

    The header lists only the sample names, while every data row starts with
    the splice-site ID; the trailing separator of the header line is dropped.
    """
    with open(input_file) as f:
        header = f.readline().strip().split()
    if header[-1] == "":
        header = header[:-1]
    df = pd.read_csv(input_file, sep=" ", skiprows=1, header=None)
    df.index = df.iloc[:, 0]          # splice-site ID
    df_data = df.iloc[:, 1:]          # PSI values
    if len(header) != df_data.shape[1]:
        raise ValueError(
            f"Number of header columns ({len(header)}) does not match number "
            f"of data columns ({df_data.shape[1]}).")
    df_data.columns = header
    return df_data


def rename_columns(df):
    """Rename meta-cell columns to '<donor>:<index within donor>'."""
    donor_counter = defaultdict(int)
    col_mapping = {}
    for col in df.columns:
        # Sample names look like <library>_<donor>_<meta-cell>; the donor ID is
        # the H-prefixed field.
        match = re.search(r'_(H\d+)_', col)
        if not match:
            raise ValueError(f"Cannot extract donor ID from column name: {col}")
        donor = match.group(1)
        donor_counter[donor] += 1
        col_mapping[col] = f"{donor}:{donor_counter[donor]}"
    df.columns = [col_mapping[col] for col in df.columns]
    return df, col_mapping


def read_reference_ids(reference_file):
    """Splice-site IDs of the pseudobulk (filtered) phenotype matrix."""
    with open(reference_file) as f:
        ref_ids = [line.strip().split()[0] for line in f.readlines()[1:]]
    return set(ref_ids)


def main(input_file, output_file, mapping_file, reference_file=None):
    df = read_and_fix_header(input_file)

    if reference_file:
        ref_ids = read_reference_ids(reference_file)
        df = df[df.index.isin(ref_ids)]

    df, col_mapping = rename_columns(df)

    # Write in the same layout as the input (header = sample names only).
    with open(output_file, "w") as f:
        f.write(" ".join(df.columns) + " \n")
        for idx, row in df.iterrows():
            row_str = " ".join(map(str, row.values))
            f.write(f"{idx} {row_str}\n")

    mapping_df = pd.DataFrame(list(col_mapping.items()),
                              columns=["OriginalName", "NewName"])
    mapping_df.to_csv(mapping_file, sep=" ", index=False)


def celltype_list(ref_dir, depot):
    """Cell types tested in one depot (first column of the reference list)."""
    path = os.path.join(ref_dir, f"celltypes.QTL.{depot}.txt")
    if not os.path.exists(path):
        raise SystemExit(f"Cell-type list not found: {path}")
    with open(path) as f:
        return [line.split("\t")[0].strip() for line in f if line.strip()]


def run_batch(cell_type_filter=None, depot_filter=None):
    """Process every cell type and depot (paths relative to this script)."""
    here = os.path.dirname(os.path.abspath(__file__))
    pheno_dir = os.path.join(here, "..", "data", "04_phenotype")
    ref_dir = os.path.join(here, "..", "reference")

    depots = [depot_filter] if depot_filter else ["SAT", "VAT"]

    for depot in depots:
        cell_types = ([cell_type_filter] if cell_type_filter
                      else celltype_list(ref_dir, depot))

        for cell_type in cell_types:
            d = os.path.join(pheno_dir, "phenotype_output", cell_type, depot)
            if not os.path.isdir(d):
                print(f"skip {cell_type} ({depot}): {d} not found")
                continue

            print(f"=== {cell_type} ({depot}) ===")

            # PSI matrix
            main(os.path.join(d, "phenotype_output.noFilter.prop"),
                 os.path.join(d, "phenotype_output.noFilter_rename.prop"),
                 os.path.join(d, "mapping.noFilter.txt"),
                 os.path.join(d, "phenotype_output.all.prop"))

            # Model input matrix
            main(os.path.join(d, "phenotype_output.noFilter"),
                 os.path.join(d, "phenotype_output.noFilter_rename.filtered"),
                 os.path.join(d, "mapping2.noFilter.txt"),
                 os.path.join(d, "phenotype_output.all.filtered"))


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Filter splice sites and rename meta-cell columns")
    parser.add_argument("-i", "--input", help="Input phenotype matrix")
    parser.add_argument("-o", "--output", help="Output phenotype matrix")
    parser.add_argument("-m", "--mapping",
                        help="Output file with the old -> new column names")
    parser.add_argument("-r", "--reference",
                        help="Optional reference matrix used to filter rows")
    parser.add_argument("--cell_type", help="Restrict the batch run to one cell type")
    parser.add_argument("--depot", help="Restrict the batch run to one depot")
    args = parser.parse_args()

    if args.input:
        if not (args.output and args.mapping):
            parser.error("-i also requires -o and -m")
        main(args.input, args.output, args.mapping, args.reference)
    else:
        run_batch(args.cell_type, args.depot)
