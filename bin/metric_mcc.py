#!/usr/bin/env python

"""
Run the Matthews Correlation Coefficient (MCC) metric

Usage:
    metric_F1.py --out-file=<path> --dataset=<str> --method=<str> [options] <file>

Options:
    -h --help            Show this screen.
    --out-file=<path>    Path to output file.
    --dataset=<str>      Name of the dataset.
    --method=<str>       Name of the method.
    --labels=<str>       Column of obs containing cell labels [default: Label].
"""

import anndata as ad
import pandas as pd
from sklearn.metrics import matthews_corrcoef

def run_MCC(adata, dataset, labels, method):

    print("Calculating MCC...")
    score = matthews_corrcoef(adata.obs[labels], adata.obs["ClusterMatched"])

    results = pd.DataFrame(
        {
            "Dataset": [dataset],
            "Method": [method],
            "Metric": ["MCC"],
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
    results = run_MCC(adata, dataset, labels, method)
    print(f"Writing data to '{out_file}'...")
    results.to_csv(out_file, sep="\t", index=False)
