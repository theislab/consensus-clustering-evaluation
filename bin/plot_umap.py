#!/usr/bin/env python

"""
Plot UMAP

Usage:
    plot_umap.py --out-file=<path> --clusters=<path> [options] <file>

Options:
    -h --help            Show this screen.
    --out-file=<path>    Path to output file.
    --clusters=<path>    Path to file containing clusters.
"""

import scanpy as sc
import pandas as pd

def plot_umap(adata, clusters):

    print("Storing clusters...")
    clusters = clusters.set_index("Cell")
    adata.obs["Cluster"] = clusters["Cluster"].astype("category")

    print("Creating UMAP plot...")
    umap = sc.pl.umap(
        adata,
        color=["Label", "Cluster"],
        legend_fontsize="small",
        legend_fontweight="light",
        add_outline=True,
        outline_width=(0.1, 0.05),
        ncols=1,
        show=False,
        return_fig=True
    )

    return umap

if __name__=="__main__":
    from docopt import docopt

    args = docopt(__doc__)

    file = args["<file>"]
    out_file = args["--out-file"]
    clusters_file = args["--clusters"]

    print(f"Reading data from '{file}'...")
    adata = sc.read_h5ad(file)
    print("Read data:")
    print(adata)
    print(f"Reading clusters from '{file}'...")
    clusters = pd.read_csv(clusters_file, sep="\t")
    print("Read clusters:")
    print(clusters)
    umap = plot_umap(adata, clusters)
    print(f"Writing plot to '{out_file}'...")
    umap.savefig(out_file, dpi=300, bbox_inches="tight")
    print("Done!")
