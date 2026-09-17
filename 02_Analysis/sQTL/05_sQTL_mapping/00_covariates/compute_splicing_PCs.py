#!/usr/bin/env python3
"""PCA of a splicing phenotype matrix.

Compute the top principal components of the PSI matrix of one cell type and
depot. They are used as covariates in the cis-sQTL model to absorb unwanted
variation between meta-cells (technical and unmeasured biological effects).
"""

import argparse

import matplotlib.pyplot as plt
import pandas as pd
from sklearn.decomposition import PCA


def run_pca(input_file, output_file, n_components=10, plot_file=None):
    # Splice sites in rows, samples in columns
    df = pd.read_csv(input_file, sep=r"\s+", index_col=0)

    # PCA over the samples, so the matrix is transposed first
    X = df.T
    pca = PCA(n_components=n_components)
    principal_components = pca.fit_transform(X)

    pc_df = pd.DataFrame(
        principal_components,
        columns=[f'PC{i + 1}' for i in range(principal_components.shape[1])],
        index=X.index)

    pc_df.to_csv(output_file, sep="\t")

    if plot_file:
        plt.figure(figsize=(8, 6))
        plt.scatter(pc_df['PC1'], pc_df['PC2'])
        plt.xlabel("PC1")
        plt.ylabel("PC2")
        plt.title("PCA of Samples")
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(plot_file)
        plt.close()

    print("PCA complete. Results saved to:", output_file)
    if plot_file:
        print("Plot saved to:", plot_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run PCA on a splicing phenotype matrix.")
    parser.add_argument("--input", required=True,
                        help="PSI matrix (splice sites x samples)")
    parser.add_argument("--output", required=True,
                        help="Output table of principal components")
    parser.add_argument("--n_components", type=int, default=5,
                        help="Number of principal components to compute")
    parser.add_argument("--plot",
                        help="Optional path for a PC1 vs PC2 scatter plot")
    args = parser.parse_args()

    run_pca(args.input, args.output, args.n_components, args.plot)
