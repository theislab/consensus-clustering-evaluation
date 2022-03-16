#!/usr/bin/env python

"""
Prepare dataset

Usage:
    prep_dataset.py --name=<str> --labels=<str> [options] <file>

Options:
    -h --help        Show this screen.
    --name=<str>     Name of the dataset.
    --labels=<str>   Column of obs containing cell labels.
"""

import scanpy as sc

def prep_dataset(name, file, labels):
    adata = sc.read_h5ad(file)

    print(name)
    print(f"File: {file}")
    print(f"Cells: {adata.n_obs}")
    print(f"Genes: {adata.n_vars}")
    print(f"Labels: {labels}")


if __name__=="__main__":
    from docopt import docopt

    args = docopt(__doc__)

    file = args["<file>"]
    name = args["--name"]
    labels = args["--labels"]

    adata = prep_dataset(name, file, labels)
