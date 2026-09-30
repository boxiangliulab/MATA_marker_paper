#!/usr/bin/env python3
import argparse
import pandas as pd
from pyliftover import LiftOver

def main():
    ap = argparse.ArgumentParser(description="LiftOver hg19 table (chr,pos,ref,alt,SNP,...) to hg38 using UCSC chain.")
    ap.add_argument("--input", required=True, help="Input TSV (hg19), must include columns: chr, pos, ref, alt, SNP")
    ap.add_argument("--outdir", required=True, help="Output directory")
    ap.add_argument("--chain", required=True, help="hg19ToHg38.over.chain.gz (UCSC chain file)")
    ap.add_argument("--outfile", default="All_chr1_22_merged.hg38.txt", help="Output file name")
    args = ap.parse_args()
    # Read table
    df = pd.read_csv(args.input, sep="\t", dtype={"chr": str})
    required = {"chr", "pos", "ref", "alt", "SNP"}
    missing = required - set(df.columns)
    if missing:
        raise SystemExit(f"[ERROR] Missing required columns: {sorted(missing)}")
    df["pos"] = df["pos"].astype(int)
    # Liftover
    lo = LiftOver(args.chain)
    chr38_list = []
    pos38_list = []
    mapped_list = []
    for c, p in zip(df["chr"].tolist(), df["pos"].tolist()):
        chrom = c
        if not chrom.startswith("chr"):
            chrom = "chr" + chrom
        res = lo.convert_coordinate(chrom, p - 1)
        if not res:
            chr38_list.append(None)
            pos38_list.append(None)
            mapped_list.append(False)
            continue
        # Take the first mapping
        new_chr, new_pos0, strand, score = res[0]
        # Convert back to 1-based
        new_pos1 = int(new_pos0) + 1
        chr38_list.append(new_chr.replace("chr", ""))
        pos38_list.append(new_pos1)
        mapped_list.append(True)
    df["chr38"] = chr38_list
    df["pos38"] = pos38_list
    df["mapped_hg38"] = mapped_list
    unmapped = df.loc[~df["mapped_hg38"], ["chr", "pos", "ref", "alt", "SNP"]].copy()
    mapped = df.loc[df["mapped_hg38"]].copy()
    # Create hg38 variant id in GTEx style
    mapped["SNP38"] = (
        mapped["chr38"].astype(str) + "_" +
        mapped["pos38"].astype(int).astype(str) + "_" +
        mapped["ref"].astype(str) + "_" +
        mapped["alt"].astype(str)
    )
    # Write outputs
    out_main = f"{args.outdir.rstrip('/')}/{args.outfile}"
    out_unmapped = f"{args.outdir.rstrip('/')}/All_chr1_22_merged.hg38.unmapped.tsv"
    mapped.to_csv(out_main, sep="\t", index=False)
    unmapped.to_csv(out_unmapped, sep="\t", index=False)
    print(f"[OK] Mapped: {mapped.shape[0]} / {df.shape[0]}")
    print(f"[OK] Output: {out_main}")
    print(f"[OK] Unmapped report: {out_unmapped}")

if __name__ == "__main__":
    main()
