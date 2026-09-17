#!/usr/bin/env python3
"""Turn the meta-cell map into per-sample barcode lists.

The junction files produced in stage 02 are indexed by library, so the
meta-cell assignment (a single CSV per cell type and depot) has to be
broken up into one barcode list per library / donor. Two groupings are
produced:

  --mode meta_cell : one file per meta-cell, used for the meta-cell level
                     (single-cell resolution) sQTL phenotypes
                         <out_dir>/<depot>/<library>_<donor>/<meta_cell>.txt
  --mode pseudobulk: one file per donor, i.e. all nuclei of that donor and
                     depot pooled into one pseudobulk sample, used for the
                     pseudobulk statistics that filter the splice sites
                         <out_dir>/all/<depot>/<library>_<donor>/<donor>.txt

Barcodes in the meta-cell map carry a library suffix (<barcode>_B<n>_L<n>).
The suffix identifies the library and is stripped from the barcode itself,
because the junction files store the plain barcode.
"""

import argparse
import csv
import os
import re
from collections import defaultdict

parser = argparse.ArgumentParser(
    description="Split a meta-cell map into per-sample barcode lists")
parser.add_argument("--input", nargs="+", required=True,
                    help="Meta-cell CSV file(s) from 02_assign_meta_cells_knn.py")
parser.add_argument("--out_dir", required=True,
                    help="Output directory (one sub-directory per depot)")
parser.add_argument("--mode", choices=["meta_cell", "pseudobulk"],
                    required=True,
                    help="Group barcodes by meta-cell or into one pseudobulk\n"
                         "sample per donor")
args = parser.parse_args()

# Library suffix embedded in the barcode, e.g. AACCTG..._B12_L3
LIBRARY_PATTERN = re.compile(r"_B\d+_L\d+")

barcodes_by_file = defaultdict(list)

for infile in args.input:
    with open(infile) as f:
        for row in csv.DictReader(f):
            raw_cell = row["cell"]

            match = LIBRARY_PATTERN.search(raw_cell)
            if not match:
                continue                       # skip malformed entries
            library = match.group(0)[1:]       # strip the leading underscore
            barcode = LIBRARY_PATTERN.sub("", raw_cell)

            donor = row["donor"]
            depot = row["depot"]
            sample = f"{library}_{donor}"      # matches the junction file names

            if args.mode == "meta_cell":
                out_path = os.path.join(
                    args.out_dir, depot, sample, f"{row['meta_cell']}.txt")
            else:
                out_path = os.path.join(
                    args.out_dir, "all", depot, sample, f"{donor}.txt")

            barcodes_by_file[out_path].append(barcode)

for path, barcodes in barcodes_by_file.items():
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        f.write("\n".join(barcodes) + "\n")

print(f"Wrote {len(barcodes_by_file)} barcode lists under {args.out_dir}")
