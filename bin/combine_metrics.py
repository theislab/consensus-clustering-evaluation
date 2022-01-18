#!/usr/bin/env python

"""
Combine metrics output files

Usage:
    combine_metrics.py --out-file=<path> [options] <file>...

Options:
    -h --help            Show this screen.
    --out-file=<path>    Path to output file.
"""

import pandas as pd

def combine_metrics(files):

    dfs = (pd.read_csv(file, sep="\t") for file in files)
    metrics = pd.concat(dfs, ignore_index=True)

    return metrics

if __name__=="__main__":
    from docopt import docopt

    args = docopt(__doc__)

    files = args["<file>"]
    out_file = args["--out-file"]

    metrics = combine_metrics(files)

    print(f"Writing data to '{out_file}'...")
    metrics.to_csv(out_file, sep="\t", index=False)
