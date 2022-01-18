#!/usr/bin/env python

"""
Run the scanpy method

Usage:
    method_scanpy.py --out-file=<path> [options] <file>

Options:
    -h --help            Show this screen.
    --out-file=<path>    Path to output file.
"""

import scanpy as sc

def run_scanpy(adata):

    print("Scaling counts to counts per 10000...")
    sc.pp.normalize_total(adata, target_sum=1e4)
    print("Log transforming...")
    sc.pp.log1p(adata)
    print("Identifying highly variable genes...")
    sc.pp.highly_variable_genes(adata)
    print(f"Found {adata.var['highly_variable'].sum()} HVGs")
    print("Calculating PCA...")
    sc.tl.pca(adata, svd_solver="arpack")
    print(f"Building KNN graph...")
    sc.pp.neighbors(adata)
    print(f"Running Leiden community detection...")
    sc.tl.leiden(adata)
    adata.obs["Cluster"] = adata.obs["leiden"]

    return adata

if __name__=="__main__":
    from docopt import docopt

    args = docopt(__doc__)

    file = args["<file>"]
    out_file = args["--out-file"]

    print(f"Reading data from '{file}'...")
    adata = sc.read_h5ad(file)
    print("Read data:")
    print(adata)
    adata = run_scanpy(adata)
    print(f"Writing data to '{out_file}'...")
    adata.write_h5ad(out_file)
