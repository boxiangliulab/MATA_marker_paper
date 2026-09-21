#!/usr/bin/env python3
"""Extract VAT CellPhoneDB interaction summaries.

Inputs are provided through SECTION4_VAT_H5AD, SECTION4_CELLPHONEDB_H5AD, and
SECTION4_CELLPHONEDB_PKL. Outputs include sender-receiver count matrices,
interaction-weight matrices, and edge tables with pathway-strength summaries.
The exported tables include DP-ATM IL1B-receptor interactions used in the
downstream analysis.
"""
import os
import numpy as np
import pandas as pd
import omicverse as ov
import pickle


# =====================================================
# 0) paths (VAT only)
# =====================================================
adata_path = os.environ["SECTION4_VAT_H5AD"]
cpdb_h5ad = os.environ["SECTION4_CELLPHONEDB_H5AD"]
cpdb_pkl = os.environ["SECTION4_CELLPHONEDB_PKL"]
out_root = os.path.join(os.environ["SECTION4_OUTPUT"], "cellphonedb_interactions")
output_dir = os.path.join(out_root, "tsv")
os.makedirs(output_dir, exist_ok=True)
print("[info] output directory:", output_dir)

# =====================================================
# 1) load h5ad
# =====================================================
adata = ov.read(adata_path)
adata_cpdb = ov.read(cpdb_h5ad)

# Standardize sender and receiver labels for downstream table names.
for col in ["sender", "receiver"]:
    if col in adata_cpdb.obs.columns:
        adata_cpdb.obs[col] = adata_cpdb.obs[col].astype(str).str.replace("+", "plus", regex=False)
    else:
        raise RuntimeError(f"adata_cpdb.obs missing column '{col}'")

# Initialize the result helper used to calculate interaction summaries.
interaction_helper = ov.pl.CellChatViz(adata_cpdb, palette=None)

# =====================================================
# 2) thresholds
# =====================================================
pvalue_threshold = 0.05
mean_threshold   = 0.1   # Edge summaries require p <= 0.05 and mean >= 0.1.
min_expression   = 0.1   # Pathway strength uses mean >= 0.1 without a p-value filter.
pathway_method   = "mean"

target_pathways = [
    "Signaling by Interleukin",
    "Signaling by Plasminogen Activator",
    "Adhesion by ICAM",
]

def safe_name(s: str) -> str:
    return (
        s.replace("/", "_")
         .replace(" ", "_")
         .replace(":", "_")
         .replace("-", "_")
    )

# =====================================================
# 3) means/pvalues (from h5ad layers)
# =====================================================
def _to_np(x):
    return x.A if hasattr(x, "A") else np.asarray(x)

for k in ["means", "pvalues"]:
    if k not in adata_cpdb.layers:
        raise RuntimeError(f"adata_cpdb.layers missing '{k}'. layers={list(adata_cpdb.layers.keys())}")

means_np = _to_np(adata_cpdb.layers["means"])
pvals_np = _to_np(adata_cpdb.layers["pvalues"])

# =====================================================
# 4) interaction_scores (from pkl; LR × (meta+cellpair cols))
#    -> obs(cellpair) × var(LR)
# =====================================================
print("[info] loading pkl:", cpdb_pkl)
with open(cpdb_pkl, "rb") as f:
    pkl = pickle.load(f)

if "interaction_scores" not in pkl:
    raise RuntimeError(f"pkl missing 'interaction_scores'. keys={list(pkl.keys())}")

df_scores = pkl["interaction_scores"]
if not isinstance(df_scores, pd.DataFrame):
    raise RuntimeError(f"pkl['interaction_scores'] must be DataFrame. got={type(df_scores)}")

meta_cols = [
    "id_cp_interaction","interacting_pair","partner_a","partner_b","gene_a","gene_b",
    "secreted","receptor_a","receptor_b","annotation_strategy","is_integrin",
    "directionality","classification"
]
cellpair_cols = [c for c in df_scores.columns if c not in meta_cols]
if len(cellpair_cols) == 0:
    raise RuntimeError("No cellpair columns found in pkl interaction_scores after removing meta cols.")

# h5ad obs cellpair label: "sender|receiver"
obs_cellpairs = (adata_cpdb.obs["sender"].astype(str) + "|" + adata_cpdb.obs["receiver"].astype(str)).tolist()

# h5ad var LR label
if "interacting_pair" in adata_cpdb.var.columns:
    var_lr = adata_cpdb.var["interacting_pair"].astype(str).tolist()
else:
    var_lr = adata_cpdb.var_names.astype(str).tolist()

pkl_lr = df_scores["interacting_pair"].astype(str).tolist()

# 1) cellpair columns align
missing_cp = sorted(set(obs_cellpairs) - set(cellpair_cols))
if len(missing_cp) > 0:
    raise RuntimeError(
        "pkl interaction_scores is missing some cellpair columns.\n"
        f"Example missing: {missing_cp[:10]}\n"
        "Check sender/receiver naming ('+' vs 'plus')."
    )

df_scores_cp = df_scores.loc[:, cellpair_cols]  # LR × cellpair

# 2) LR rows align by name if possible
if set(var_lr).issubset(set(pkl_lr)):
    lr2idx = {}
    for i, lr in enumerate(pkl_lr):
        if lr not in lr2idx:  # first occurrence only
            lr2idx[lr] = i
    row_idx = [lr2idx[lr] for lr in var_lr]
    df_scores_cp = df_scores_cp.iloc[row_idx, :]
else:
    if df_scores_cp.shape[0] != len(var_lr):
        raise RuntimeError(
            f"LR dimension mismatch: pkl LR={df_scores_cp.shape[0]} vs h5ad var={len(var_lr)} "
            "and cannot align by name."
        )
    print("[warn] could not align LR rows by name. Using raw row order (assumes same order).")

# 3) reorder columns to obs order, then transpose -> obs×var
df_scores_cp = df_scores_cp.loc[:, obs_cellpairs]   # LR × cellpair (obs order)
scores_np = df_scores_cp.to_numpy().T               # cellpair × LR

print("[info] scores_np shape:", scores_np.shape, "expected:", means_np.shape)
if scores_np.shape != means_np.shape:
    raise RuntimeError(f"scores_np shape {scores_np.shape} != means_np shape {means_np.shape} (must match obs×var).")

# =====================================================
# 5) cell types mapping
# =====================================================
cell_types = interaction_helper._get_unique_cell_types()
ct2i = {ct: i for i, ct in enumerate(cell_types)}
senders   = adata_cpdb.obs["sender"].astype(str).tolist()
receivers = adata_cpdb.obs["receiver"].astype(str).tolist()

# =====================================================
# 6) Calculate pathway summaries.
# =====================================================
def _agg_strength(x: np.ndarray, method: str):
    if x.size == 0:
        return 0.0
    if method == "mean":
        return float(np.mean(x))
    if method == "sum":
        return float(np.sum(x))
    if method == "max":
        return float(np.max(x))
    if method == "median":
        return float(np.median(x))
    raise ValueError(f"Unknown method: {method}")

def calc_pathway_mats(pathway_name: str,
                      pthr: float,
                      mthr: float,
                      min_expr: float,
                      method: str):
    col_mask = (adata_cpdb.var["classification"].astype(str) == pathway_name).values
    cols = np.where(col_mask)[0]
    if len(cols) == 0:
        return None, None, None

    n = len(cell_types)
    count_mat = np.zeros((n, n), dtype=np.int64)
    weight_mat = np.zeros((n, n), dtype=np.float64)
    score_sum_mat = np.zeros((n, n), dtype=np.float64)

    M = means_np[:, cols]
    P = pvals_np[:, cols]
    S = scores_np[:, cols]

    # A. Significant interactions passing the mean threshold.
    keep_sig = (P <= pthr) & (M >= mthr)
    row_counts_sig = keep_sig.sum(axis=1)
    row_weight_means_sig = (M * keep_sig).sum(axis=1)
    row_score_sum_sig = (S * keep_sig).sum(axis=1)

    # B. Pathway strength as implemented by the source helper.
    valid_mask = (M >= min_expr)

    row_strength = np.zeros((adata_cpdb.n_obs,), dtype=np.float64)
    row_pathway_pval = np.ones((adata_cpdb.n_obs,), dtype=np.float64)
    row_n_valid = np.zeros((adata_cpdb.n_obs,), dtype=np.int64)

    # (C) sigOnly strength: (mean >= min_expr) & (p <= pthr)
    valid_sig_mask = (M >= min_expr) & (P <= pthr)
    row_strength_sigOnly = np.zeros((adata_cpdb.n_obs,), dtype=np.float64)

    for r in range(adata_cpdb.n_obs):
        # Default pathway-strength definition.
        vm = valid_mask[r, :]
        if np.any(vm):
            valid_means = M[r, vm]
            valid_pvals = P[r, vm]
            row_strength[r] = _agg_strength(valid_means, method=method)
            row_pathway_pval[r] = float(np.min(valid_pvals))
            row_n_valid[r] = int(valid_means.size)
        else:
            row_strength[r] = 0.0
            row_pathway_pval[r] = 1.0
            row_n_valid[r] = 0

        # Alternative strength restricted to significant interactions.
        vms = valid_sig_mask[r, :]
        if np.any(vms):
            valid_means_sig = M[r, vms]
            row_strength_sigOnly[r] = _agg_strength(valid_means_sig, method=method)
        else:
            row_strength_sigOnly[r] = 0.0

    # Accumulate matrices used in edge summaries.
    for r in range(adata_cpdb.n_obs):
        s = senders[r]
        t = receivers[r]
        if (s not in ct2i) or (t not in ct2i):
            continue
        i = ct2i[s]
        j = ct2i[t]

        c = int(row_counts_sig[r])
        if c > 0:
            count_mat[i, j] += c
            weight_mat[i, j] += float(row_weight_means_sig[r])
            score_sum_mat[i, j] += float(row_score_sum_sig[r])

    # Aggregate by sender-receiver pair.
    key_strength_vals, key_pval_vals, key_nvalid_vals, key_strength_sigOnly_vals = {}, {}, {}, {}
    for r in range(adata_cpdb.n_obs):
        s = senders[r]
        t = receivers[r]
        if (s not in ct2i) or (t not in ct2i):
            continue
        key = (s, t)
        key_strength_vals.setdefault(key, []).append(float(row_strength[r]))
        key_pval_vals.setdefault(key, []).append(float(row_pathway_pval[r]))
        key_nvalid_vals.setdefault(key, []).append(int(row_n_valid[r]))
        key_strength_sigOnly_vals.setdefault(key, []).append(float(row_strength_sigOnly[r]))

    edges = []
    for i, s in enumerate(cell_types):
        for j, t in enumerate(cell_types):
            c = int(count_mat[i, j])
            if c == 0:
                continue

            w = float(weight_mat[i, j])
            sc = float(score_sum_mat[i, j])

            key = (s, t)
            if key in key_strength_vals:
                pw_strength = float(np.mean(key_strength_vals[key]))
                pw_pval = float(np.min(key_pval_vals[key]))
                pw_nvalid = int(np.sum(key_nvalid_vals[key]))
                pw_strength_sigOnly = float(np.mean(key_strength_sigOnly_vals[key]))
            else:
                pw_strength, pw_pval, pw_nvalid, pw_strength_sigOnly = 0.0, 1.0, 0, 0.0

            edges.append({
                "pathway": pathway_name,
                "sender": s,
                "receiver": t,

                # Summary based on significance and mean thresholds.
                "count_sig_pairs": c,
                "weight_sum_means_sig": w,
                "mean_per_sig_pair": (w / c) if c > 0 else np.nan,

                "interaction_score_sum_sig": sc,
                "interaction_score_mean_per_sig_pair": (sc / c) if c > 0 else np.nan,

                # Source-method pathway strength.
                f"pathway_strength_{method}_minExpr{min_expr}": pw_strength,
                "pathway_pval_minValid": pw_pval,
                f"n_valid_interactions_minExpr{min_expr}": pw_nvalid,

                # Additional significance-restricted strength.
                f"pathway_strength_sigOnly_{method}": pw_strength_sigOnly,
            })

    return count_mat, weight_mat, pd.DataFrame(edges)


# =====================================================
# 7) run + save
# =====================================================
for pw in target_pathways:
    pw_tag = safe_name(pw)

    res = calc_pathway_mats(
        pw,
        pthr=pvalue_threshold,
        mthr=mean_threshold,
        min_expr=min_expression,
        method=pathway_method
    )
    if res[0] is None:
        print("[skip] no LR in pathway:", pw)
        continue

    count_mat_pw, weight_mat_pw, edges_df_pw = res
    count_df_pw  = pd.DataFrame(count_mat_pw, index=cell_types, columns=cell_types)
    weight_df_pw = pd.DataFrame(weight_mat_pw, index=cell_types, columns=cell_types)

    count_df_pw.to_csv(
        os.path.join(
            output_dir,
            f"VAT__pathway__{pw_tag}__count_matrix__p{pvalue_threshold}"
            f"__mean{mean_threshold}.tsv",
        ),
        sep="\t"
    )
    weight_df_pw.to_csv(
        os.path.join(
            output_dir,
            f"VAT__pathway__{pw_tag}__weight_matrix_sumMeansSig"
            f"__p{pvalue_threshold}__mean{mean_threshold}.tsv",
        ),
        sep="\t"
    )
    edges_df_pw.to_csv(
        os.path.join(
            output_dir,
            f"VAT__pathway__{pw_tag}__edges__p{pvalue_threshold}"
            f"__mean{mean_threshold}__minExpr{min_expression}"
            f"__method{pathway_method}.tsv",
        ),
        sep="\t",
        index=False
    )

    print("[saved]", pw, "| edge rows:", edges_df_pw.shape[0])

print("[done] TSV files:", output_dir)
