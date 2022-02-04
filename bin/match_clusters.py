#!/usr/bin/env python

"""
Match cluster labels using the Hungarian algorithm

Usage:
    match_clusters.py --out-file=<path> --labels=<str> --clusters=<str> [options] <file>

Options:
    -h --help            Show this screen.
    --out-file=<path>    Path to output file.
    --labels=<str>       Column of obs containing cell labels.
    --clusters=<str>     Column of obs containing cluster labels.
"""

import anndata as ad
import numpy as np
from scipy.spatial import distance
from scipy.optimize import linear_sum_assignment

def match_clusters(adata, labels_col, clusters_col):

    print(f"Matching clusters in '{clusters_col}' to labels in '{labels_col}'...")

    adata.obs[labels_col] = adata.obs[labels_col].astype("category")
    labels = list(adata.obs[labels_col])
    label_levels = list(adata.obs[labels_col].cat.categories)

    adata.obs[clusters_col] = adata.obs[clusters_col].astype("category")
    clusters = list(adata.obs[clusters_col])
    cluster_levels = list(adata.obs[clusters_col].cat.categories)

    jaccard = np.zeros((len(label_levels), len(cluster_levels)))
    combos = [(l, c) for l in label_levels for c in cluster_levels]

    for label, cluster in combos:
        labels_bin = [1 if l == label else 0 for l in labels]
        clusters_bin = [1 if c == cluster else 0 for c in clusters]

        label_idx = label_levels.index(label)
        cluster_idx = cluster_levels.index(cluster)
        jaccard[label_idx, cluster_idx] = distance.jaccard(labels_bin, clusters_bin)

    map = {}
    while not all(cluster in map for cluster in cluster_levels):
        not_matched = [c for c in cluster_levels if c not in map.keys()]
        print(f"Unmatched clusters: {', '.join(not_matched)}")
        not_matched_idx = [cluster_levels.index(c) for c in not_matched]
        matches = linear_sum_assignment(jaccard[:, not_matched_idx])

        print("Found these matches:")
        for label, cluster in zip(matches[0], matches[1]):
            cluster_level = not_matched[cluster]
            label_level = label_levels[label]
            print(f"\tCluster {cluster_level} -> Label {label_level}")
            map[cluster_level] = label_level

    print("All clusters matched")
    adata.obs["ClusterMatched"] = adata.obs[clusters_col]
    adata.obs = adata.obs.replace(dict(ClusterMatched = map))

    return adata

if __name__=="__main__":
    from docopt import docopt

    args = docopt(__doc__)

    file = args["<file>"]
    out_file = args["--out-file"]
    labels = args["--labels"]
    clusters = args["--clusters"]

    print(f"Reading data from '{file}'...")
    adata = ad.read_h5ad(file)
    print("Read data:")
    print(adata)
    adata = match_clusters(adata, labels, clusters)
    print(f"Writing data to '{out_file}'...")
    adata.write_h5ad(out_file)
