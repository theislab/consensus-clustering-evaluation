#!/usr/bin/env python

"""
Run the random method

Usage:
    method_random.py --out-file=<path> --labels=<str> [options] <file>

Options:
    -h --help            Show this screen.
    --labels=<str>       Column of obs containing cell labels.
    --out-file=<path>    Path to output file.
"""

import anndata
import numpy as np

def run_random(adata, labels):

    print("Generating random cluster levels...")
    label_levels = list(adata.obs[labels].cat.categories)

    rng = np.random.default_rng(seed=1)

    adata.obs["Cluster"] = rng.choice(len(label_levels), adata.n_obs)
    adata.obs["Cluster"] = adata.obs["Cluster"].astype('category')

    return adata

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
    adata = run_random(adata, labels)
    print(f"Writing data to '{out_file}'...")
    adata.write_h5ad(out_file)
