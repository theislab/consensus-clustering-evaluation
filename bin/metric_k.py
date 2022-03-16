#!/usr/bin/env python

"""
Run the number of clusters (k) metric

Usage:
    metric_k.py --out-file=<path> --dataset=<str> --method=<str> [options] <file>

Options:
    -h --help            Show this screen.
    --out-file=<path>    Path to output file.
    --dataset=<str>      Name of the dataset.
    --method=<str>       Name of the method.
    --labels=<str>       Column containing cell labels [default: Label].
"""

import pandas as pd

def run_k(clusters, dataset, labels, method):

    print("Calculating k...")
    score = clusters["Cluster"].nunique()
    prop = score / clusters[labels].nunique()

    results = pd.DataFrame(
        {
            "Dataset": [dataset, dataset],
            "Method": [method, method],
            "Metric": ["k", "k_prop"],
            "Score": [score, prop]
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
    clusters = pd.read_csv(file, sep="\t")
    print("Read data:")
    print(clusters)
    results = run_k(clusters, dataset, labels, method)
    print(f"Writing data to '{out_file}'...")
    results.to_csv(out_file, sep="\t", index=False)
