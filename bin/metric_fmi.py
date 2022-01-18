#!/usr/bin/env python

"""
Run the Fowlkes-Mallows Index (FMI) metric

Usage:
    metric_fmi.py --out-file=<path> --dataset=<str> --labels=<str> --method=<str> [options] <file>

Options:
    -h --help            Show this screen.
    --out-file=<path>    Path to output file.
    --dataset=<str>      Name of the dataset.
    --labels=<str>       Column of obs containing cell labels.
    --method=<str>       Name of the method.
"""

import anndata as ad
import pandas as pd
from sklearn.metrics import fowlkes_mallows_score

def run_fmi(adata, dataset, labels, method):

    print("Calculating FMI...")
    score = fowlkes_mallows_score(adata.obs[labels], adata.obs["Cluster"])

    results = pd.DataFrame(
        {
            "Dataset": [dataset],
            "Method": [method],
            "Metric": ["FMI"],
            "Score": [score]
        }
    )

    return results

if __name__=="__main__":
    from docopt import docopt

    args = docopt(__doc__)

    file = args["<file>"]
    out_file = args["--out-file"]
    dataset = args["--dataset"]
    labels = args["--labels"]
    method = args["--method"]

    print(f"Reading data from '{file}'...")
    adata = ad.read_h5ad(file)
    print("Read data:")
    print(adata)
    results = run_fmi(adata, dataset, labels, method)
    print(f"Writing data to '{out_file}'...")
    results.to_csv(out_file, sep="\t", index=False)