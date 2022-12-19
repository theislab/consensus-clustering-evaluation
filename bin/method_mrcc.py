#!/usr/bin/env python

"""
Run the Multi-resolution Consensus Clustering (MRCC) method

Usage:
    method_mrcc.py --out-file=<path> [options] <file>

Options:
    -h --help                       Show this screen.
    --out-file=<path>               Path to output file.
    --graph-type=<str>              Graph type for selecting clusters (neighbour, all) [default: neighbour].
    --community-type=<str>          Multi-resolution community detection method to use (leiden, hdbscan, louvain) [default: leiden].
    --single-resolution=<float>     Parameter for the single-resolution graph community detection method (leiden resolution parameter) [default: 1].
    --multi-resolution=<float>      Parameter for the mult-resolution graph community detection method (leiden resolution parameter) [default: 1].
"""

import scanpy as sc
import pandas as pd

# Load MRCC from submodule as described here https://stackoverflow.com/a/50395128/4384120
import importlib.util
import sys
from pathlib import Path

HERE = Path(__file__).parent
MODULE_PATH = HERE.joinpath("multires-consensus-clustering/multires_consensus_clustering/__init__.py")
MODULE_NAME = "multires_consensus_clustering"
SPEC = importlib.util.spec_from_file_location(MODULE_NAME, MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)
import multires_consensus_clustering as mrcc

def run_mrcc(adata, neighbour_based, community_type, single_resolution,
             multi_resolution):

    settings = adata.uns["constclust"]["settings"]
    clusterings = adata.uns["constclust"]["clusterings"]
    clusterings = clusterings.rename_axis("cell").reset_index()

    print("Performing multi-resolution consensus clustering...")
    print(f"Neighbour-based: {neighbour_based}")
    print(f"Community type: {community_type}")
    print(f"Single-resolution parameter: {single_resolution}")
    print(f"Multi-resolution parameter: {multi_resolution}")

    print("Calculating multi-resolution graph...")
    multires_graph = mrcc.multiresolution_graph(
        clusterings,
        settings,
        "all",
        neighbour_based=neighbour_based,
        single_resolution=single_resolution
    )

    print("Clustering multi-resolution graph...")
    multires_graph = mrcc.multires_community_detection(
        multires_graph,
        community_detection=community_type,
        clustering_data=clusterings,
        multi_resolution=multi_resolution
    )

    print("Assigning cluster labels...")
    df_clusters = mrcc.graph_to_cell_labels_df(multires_graph)
    cluster_labels = mrcc.df_cell_clusters_to_labels(
        df_clusters,
        adata=adata,
        plot_labels=False
    )

    clusters = pd.DataFrame(
        {
            "Cell": adata.obs_names,
            "Label": adata.obs["Label"],
            "Cluster": cluster_labels
        }
    )

    return clusters

if __name__=="__main__":
    from docopt import docopt

    args = docopt(__doc__)

    file = args["<file>"]
    out_file = args["--out-file"]
    graph_type = args["--graph-type"]
    community_type = args["--community-type"]
    single_resolution = float(args["--single-resolution"])
    multi_resolution = float(args["--multi-resolution"])

    print(f"Reading data from '{file}'...")
    adata = sc.read_h5ad(file)
    print("Read data:")
    print(adata)

    if graph_type == "neighbour":
        neighbour_based = True
    else:
        neighbour_based = False

    clusters = run_mrcc(
        adata,
        neighbour_based,
        community_type,
        single_resolution,
        multi_resolution
    )
    print(f"Writing clusters to '{out_file}'...")
    clusters.to_csv(out_file, sep="\t", index=False)
    print("Done!")
