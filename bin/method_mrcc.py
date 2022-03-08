#!/usr/bin/env python

"""
Run the Multi-resolution Consensus Clustering (MRCC) method

Usage:
    method_mrcc.py --out-file=<path> [options] <file>

Options:
    -h --help            Show this screen.
    --out-file=<path>    Path to output file.
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

def run_mrcc(adata):

    settings = adata.uns["constclust"]["settings"]
    clusterings = adata.uns["constclust"]["clusterings"]
    clusterings = clusterings.rename_axis("cell").reset_index()

    print("Calculating multi-resolution graph...")
    multires_graph = mrcc.multiresolution_graph(
        clusterings,
        settings,
        "all",
        neighbour_based=True
    )

    print("Clustering multi-resolution graph...")
    multires_graph = mrcc.multires_community_detection(
        multires_graph,
        community_detection="leiden",
        merge_edges_threshold=1,
        outlier_detection="probability",
        outlier_detection_threshold=0.9,
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

    print(f"Reading data from '{file}'...")
    adata = sc.read_h5ad(file)
    print("Read data:")
    print(adata)
    adata.obs["Cluster"] = run_mrcc(adata)
    print(f"Writing data to '{out_file}'...")
    adata.write_h5ad(out_file)
