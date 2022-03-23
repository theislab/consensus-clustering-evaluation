#!/usr/bin/env python

"""
Run the constclust to generate multiple clusterings

Usage:
    run_constclust.py --out-file=<path> [options] <file>

Options:
    -h --help            Show this screen.
    --out-file=<path>    Path to output file.
"""

import scanpy as sc
import constclust as cc
import pandas as pd

def run_constclust(adata):

    print("Scaling counts to counts per 10000...")
    scanpy.pp.normalize_total(adata, target_sum=1e4)

    print("Log transforming...")
    scanpy.pp.log1p(adata)

    print("Finding highly variable genes...")
    sc.pp.highly_variable_genes(adata, flavor="cell_ranger", n_top_genes=2000)

    print("Performing PCA...")
    sc.tl.pca(adata)

    neighbours = [5, 15, 30, 50]
    resolutions = [0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 1.0, 2.0, 5.0, 10.0]
    states = [0, 1, 2]
    distances = ["euclidean", "correlation", "cosine"]

    results = {}

    for distance in distances:
        print(f"Clustering using {distance} distance:")
        settings, clusterings = cc.cluster(
            adata,
            n_neighbors=neighbours,
            resolutions=resolutions,
            random_state=states,
            neighbor_kwargs={"use_rep" : "X_pca", "metric" : distance},
            n_procs=1
        )
        settings["distance"] = distance
        settings["n_clusters"] = [clusterings[clustering].nunique() for clustering in clusterings.columns]
        results[distance] = {"settings" : settings, "clusterings" : clusterings}

    settings = pd.concat([results[distance]["settings"] for distance in distances])
    settings["id"] = [f"C{row:03}" for row in range(1, len(settings) + 1)]
    settings = settings[["id", "distance", "n_neighbors", "resolution", "random_state", "n_clusters"]]
    settings

    clusterings = pd.concat([results[distance]["clusterings"] for distance in distances], axis=1)
    clusterings = clusterings.set_axis(settings["id"], axis=1)
    clusterings

    adata.uns["constclust"] = {"settings" : settings, "clusterings" : clusterings}

    return adata

if __name__=="__main__":
    from docopt import docopt

    args = docopt(__doc__)

    file = args["<file>"]
    out_file = args["--out-file"]

    print(f"Reading data from '{file}'...")
    adata = sc.read_h5ad(file)
    print("Read data:")
    print(adata)
    adata = run_constclust(adata)
    print(f"Writing data to '{out_file}'...")
    adata.write_h5ad(out_file)
