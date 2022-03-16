#!/usr/bin/env python

"""
Run the random method

Usage:
    method_random.py --out-file=<path> [options] <file>

Options:
    -h --help            Show this screen.
    --labels=<str>       Column of obs containing cell labels [default: Label].
    --out-file=<path>    Path to output file.
"""

import anndata
import numpy as np
import pandas as pd

def run_random(adata, labels):

    print("Generating random cluster levels...")
    label_levels = list(adata.obs[labels].cat.categories)

    rng = np.random.default_rng(seed=1)

    cluster_labels = rng.choice(len(label_levels), adata.n_obs)

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
    labels = args["--labels"]
    out_file = args["--out-file"]

    print(f"Reading data from '{file}'...")
    adata = anndata.read_h5ad(file)
    print("Read data:")
    print(adata)
    clusters = run_random(adata, labels)
    print(f"Writing clusters to '{out_file}'...")
    clusters.to_csv(out_file, sep="\t", index=False)
    print("Done!")
