#!/usr/bin/env python3
"""Replication of the sQTL in GTEx (Storey's pi1).

Ask how many of the sQTL discovered here are also signals in the matching
GTEx tissue. Every significant (splice site, variant) pair is looked up in
the GTEx v10 cis-sQTL all-pairs files; pi1 = 1 - pi0 of the distribution of
the GTEx nominal p-values of the matched pairs is reported as the
replication rate.

Matching is done on chromosome:position of the variant and on the splice
site: GTEx variant IDs and intron coordinates are harmonised to the splice
site definition used here (chr:strand:position), with each GTEx intron
decomposed into its donor and acceptor splice site. Where several GTEx
associations match one discovery sQTL, the most significant nominal
p-value is kept. Cell types with fewer than 10 matched sQTL are excluded
from the pi1 estimate for lack of power.

Tissue: SAT -> Adipose_Subcutaneous, VAT -> Adipose_Visceral_Omentum;
it is taken from the ".SAT." / ".VAT." field of the input file name.
"""

import argparse
import os

import numpy as np
import pandas as pd

# Storey's lambda: p-values above this are assumed to come from the null.
LAMBDA = 0.5
# Below this number of matched pairs pi1 is not estimated.
MIN_MATCHED = 10


def storey_pi0(pvals, lambda_=LAMBDA):
    """Estimate the proportion of true null hypotheses (Storey's method)."""
    pvals = np.array(pvals)
    prop = np.mean(pvals > lambda_)
    pi0 = prop / (1 - lambda_)
    return min(pi0, 1.0)


def run_one(result_file, gtex_base, out_summary, out_matched):
    """Match one cell type / depot against GTEx and estimate pi1."""

    # --- Load this study's significant sQTL ----------------------------------
    top_df = pd.read_csv(result_file, sep="\t")
    top_df['snp_pos'] = top_df['SNP'].str.extract(r'^(chr[0-9XYM]+:[0-9]+)')

    # --- GTEx tissue matching the depot ------------------------------------------
    if ".SAT." in result_file:
        gtex_prefix = "Adipose_Subcutaneous"
    elif ".VAT." in result_file:
        gtex_prefix = "Adipose_Visceral_Omentum"
    else:
        raise ValueError("Unknown tissue from filename")

    matched_list = []

    # --- Match chromosome by chromosome (the GTEx files are large) ---------------
    for chr_num in range(1, 23):
        print(f"Processing chr{chr_num}...")

        fname = os.path.join(
            gtex_base, f"{gtex_prefix}.v10.cis_sqtl.allpairs.chr{chr_num}.parquet")
        if not os.path.exists(fname):
            print(f"File not found: {fname}")
            continue

        df = pd.read_parquet(fname)

        # GTEx variant IDs are chr_pos_ref_alt_b38 -> chr:pos:ref:alt
        df['SNP'] = (df['variant_id'].str.replace('_', ':', regex=False)
                                     .str.replace(':b38', '', regex=False))
        df['snp_pos'] = df['SNP'].str.extract(r'^(chr[0-9XYM]+:[0-9]+)')

        # GTEx phenotype IDs are chr:start:end:clu_<n>_<strand>; both ends of the
        # intron are splice sites and are tested separately here.
        pheno_split = df['phenotype_id'].str.split(':', expand=True)
        chrom = pheno_split[0]
        donor_pos = pheno_split[1]
        acceptor_pos = pheno_split[2]
        strand = pheno_split[3].str.split('_').str[-1]

        donor = chrom + ":" + strand + ":" + donor_pos
        acceptor = chrom + ":" + strand + ":" + acceptor_pos

        df_donor = df.copy()
        df_acceptor = df.copy()
        df_donor['splice_site'] = donor
        df_acceptor['splice_site'] = acceptor

        combined = pd.concat([df_donor, df_acceptor], ignore_index=True)
        combined.dropna(subset=['snp_pos', 'splice_site'], inplace=True)

        # Keep the most significant p-value per (splice site, variant)
        combined.sort_values(['splice_site', 'snp_pos', 'pval_nominal'],
                             inplace=True)
        combined = combined.drop_duplicates(subset=['splice_site', 'snp_pos'],
                                            keep='first')

        top_df_chr = top_df[top_df['snp_pos'].str.startswith(
            f'chr{chr_num}:')].copy()
        if top_df_chr.empty:
            continue

        combined['join_key'] = combined['splice_site'] + "|" + combined['snp_pos']
        top_df_chr['join_key'] = (top_df_chr['splice_site'] + "|"
                                  + top_df_chr['snp_pos'])

        matched_chr = pd.merge(top_df_chr, combined, on='join_key',
                               suffixes=('.result', '.gtex'))
        print(matched_chr.head())
        matched_list.append(matched_chr)

    # --- Combine and estimate pi1 -------------------------------------------------
    if matched_list:
        matched = pd.concat(matched_list, ignore_index=True)
    else:
        matched = pd.DataFrame(
            columns=['splice_site', 'SNP', 'pvalue', 'beta', 'se', 'FDR'])

    if len(matched) > MIN_MATCHED:
        pi0 = storey_pi0(matched['pval_nominal'])
        pi1 = 1 - pi0
    else:
        pi0 = np.nan
        pi1 = np.nan

    summary_df = pd.DataFrame({
        'tissue': [gtex_prefix],
        'cell_type': [os.path.basename(result_file).split('.')[0]],
        'n_matched': [len(matched)],
        'pi0': [pi0],
        'pi1': [pi1],
    })
    summary_df.to_csv(out_summary, sep='\t', index=False)

    if not matched.empty:
        matched.to_csv(out_matched, sep='\t', index=False)

    return summary_df


def celltype_list(ref_dir, depot):
    """Cell types tested in one depot (first column of the reference list)."""
    path = os.path.join(ref_dir, f"celltypes.QTL.{depot}.txt")
    if not os.path.exists(path):
        raise SystemExit(f"Cell-type list not found: {path}")
    with open(path) as f:
        return [line.split("\t")[0].strip() for line in f if line.strip()]


def main():
    parser = argparse.ArgumentParser(
        description="GTEx replication (pi1) of the sQTL")
    # Explicit single run
    parser.add_argument("--result_file", help="Significant sQTL of one cell type")
    parser.add_argument("--gtex_dir", help="Directory with the GTEx parquet files")
    parser.add_argument("--out_summary", help="Output summary table")
    parser.add_argument("--out_matched", help="Output matched associations")
    # Batch run
    parser.add_argument("--cell_type", help="Restrict the batch run to one cell type")
    parser.add_argument("--depot", help="Restrict the batch run to one depot")
    args = parser.parse_args()

    if args.result_file:
        if not all([args.gtex_dir, args.out_summary, args.out_matched]):
            parser.error("--result_file also requires --gtex_dir, --out_summary "
                         "and --out_matched")
        run_one(args.result_file, args.gtex_dir, args.out_summary,
                args.out_matched)
        return

    # --- Batch mode: every cell type and depot (paths relative to this script)
    here = os.path.dirname(os.path.abspath(__file__))
    model_dir = os.path.join(here, "..", "data", "05_sQTL_mapping")
    gtex_dir = os.path.join(here, "..", "data", "reference", "GTEx_v10_sQTL")
    out_dir = os.path.join(here, "..", "data", "06_replication")
    ref_dir = os.path.join(here, "..", "reference")
    pc_splicing = 5                      # splicing PCs, part of the file names
    g_geno = {"SAT": 10, "VAT": 6}       # genotype PCs per depot

    os.makedirs(out_dir, exist_ok=True)
    depots = [args.depot] if args.depot else ["SAT", "VAT"]
    summaries = []

    for depot in depots:
        cell_types = ([args.cell_type] if args.cell_type
                      else celltype_list(ref_dir, depot))

        for cell_type in cell_types:
            result_file = os.path.join(
                model_dir, "results",
                f"{cell_type}.{depot}.p{pc_splicing}_g{g_geno[depot]}"
                ".result.sQTL.fdr05")

            if not os.path.exists(result_file):
                print(f"skip {cell_type} ({depot}): no results")
                continue

            print(f"=== {cell_type} ({depot}) ===")
            summaries.append(run_one(
                result_file, gtex_dir,
                os.path.join(out_dir, f"{cell_type}.{depot}.pi1_summary.tsv"),
                os.path.join(out_dir, f"{cell_type}.{depot}.matched.tsv")))

    if summaries and not (args.cell_type or args.depot):
        summary_all = os.path.join(out_dir, "pi1_summary_all.tsv")
        pd.concat(summaries, ignore_index=True).to_csv(
            summary_all, sep="\t", index=False)
        print(f"Combined summary written to {summary_all}")

    print("Done.")


if __name__ == "__main__":
    main()
