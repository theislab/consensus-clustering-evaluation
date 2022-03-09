#!/usr/bin/env python

"""
Run the Multi-resolution Consensus Clustering (MRCC) method

Usage:
    method_mrcc.py --out-file=<path> [options] <file>

Options:
    -h --help                  Show this screen.
    --out-file=<path>          Path to output file.
    --neighbour-based          Whether to use neighbour-based graph connection.
    --community-type=<str>     Multi-resolution community detection method to use (component, leiden, hdbscan, louvain) [default: leiden].
    --outlier-type=<str>       Outlier detection method to use (probability, hdbscan) [default: probability].
    --outlier-thresh=<float>   Outlier detection threshold [default: 0.9].
    --merge-thresh=<float>     Edge merging threshold [default: 0.8].
"""

import scanpy as sc

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

def run_mrcc(adata, neighbour_based, community_type, outlier_type,
             outlier_thresh, merge_thresh):

    settings = adata.uns["constclust"]["settings"]
    clusterings = adata.uns["constclust"]["clusterings"]
    clusterings = clusterings.rename_axis("cell").reset_index()

    print("Performing multi-resolution consensus clustering...")
    print(f"Neighbour-based: {neighbour_based}")
    print(f"Community type: {community_type}")
    print(f"Outlier type: {outlier_type}")
    print(f"Outlier threshold: {outlier_thresh}")
    print(f"Merge threshold: {merge_thresh}")

    print("Calculating multi-resolution graph...")
    multires_graph = mrcc.multiresolution_graph(
        clusterings,
        settings,
        "all",
        neighbour_based=neighbour_based
    )

    print("Clustering multi-resolution graph...")
    multires_graph = mrcc.multires_community_detection(
        multires_graph,
        community_detection=community_type,
        merge_edges_threshold=merge_thresh,
        outlier_detection=outlier_type,
        outlier_detection_threshold=outlier_thresh,
        clustering_data=clusterings
    )

    print("Assigning cluster labels...")
    df_clusters = mrcc.graph_to_cell_labels_df(multires_graph)
    cluster_labels = mrcc.df_cell_clusters_to_labels(
        df_clusters,
        adata=adata,
        plot_labels=False
    )

    return cluster_labels

if __name__=="__main__":
    from docopt import docopt

    args = docopt(__doc__)

    file = args["<file>"]
    out_file = args["--out-file"]
    neighbour_based = args["--neighbour-based"]
    community_type = args["--community-type"]
    outlier_type = args["--outlier-type"]
    outlier_thresh = float(args["--outlier-thresh"])
    merge_thresh = float(args["--merge-thresh"])

    print(f"Reading data from '{file}'...")
    adata = sc.read_h5ad(file)
    print("Read data:")
    print(adata)
    adata.obs["Cluster"] = run_mrcc(
        adata,
        neighbour_based,
        community_type,
        outlier_type,
        outlier_thresh,
        merge_thresh
    )
    print(f"Writing data to '{out_file}'...")
    adata.write_h5ad(out_file)
