#!/usr/bin/env python

"""
Prepare dataset

Usage:
    prep_dataset.py --name=<str> --labels=<str> --out-file=<path> [options] <file>

Options:
    -h --help            Show this screen.
    --name=<str>         Name of the dataset.
    --labels=<str>       Column of obs containing cell labels.
    --out-file=<path>    Path to output file.
"""

import scanpy
import anndata

def prep_dataset(name, file, labels):
    print(f"Reading data from '{file}'...")
    adata_raw = scanpy.read_h5ad(file)

    print()
    print("====== RAW DATA =======")
    print(f"Name: {name}")
    print(f"Cells: {adata_raw.n_obs}")
    print(f"Genes: {adata_raw.n_vars}")
    print(f"Labels: {labels} ({len(adata_raw.obs[labels].cat.categories)})")
    print(adata_raw)
    print("=======================")
    print()

    print("Creating new AnnData...")
    adata = anndata.AnnData(X = adata_raw.X.copy())
    adata.obs_names = adata_raw.obs_names.copy()
    adata.var_names = adata_raw.var_names.copy()
    adata.obs["Label"] = adata_raw.obs[labels]
    del adata_raw

    print("Removing cells with less than 100 counts...")
    scanpy.pp.filter_cells(adata, min_counts=1)

    print("Removing cells with less than 100 genes...")
    scanpy.pp.filter_cells(adata, min_genes=1)

    print("Removing genes with 0 counts...")
    scanpy.pp.filter_genes(adata, min_counts=1)

    print("Removing unused labels...")
    adata.obs["Label"].cat = adata.obs["Label"].cat.remove_unused_categories()

    print("Saving counts...")
    counts = adata.X.copy()

    print("Scaling counts to counts per 10000...")
    scanpy.pp.normalize_total(adata, target_sum=1e4)

    print("Log transforming...")
    scanpy.pp.log1p(adata)

    print("Calculating 2000 highly variable genes...")
    scanpy.pp.highly_variable_genes(adata, flavor="cell_ranger", n_top_genes=2000)

    print("Calculating PCA...")
    scanpy.tl.pca(adata)

    print("Calculating neighbourhood graph...")
    scanpy.pp.neighbors(adata)

    print("Calculating UMAP...")
    scanpy.tl.umap(adata)

    print("Storing counts...")
    adata.X = counts

    print("\n==== PREPARED DATA ====")
    print(f"Name: {name}")
    print(f"Cells: {adata.n_obs}")
    print(f"Genes: {adata.n_vars}")
    print(f"Labels: Label ({len(adata.obs['Label'].cat.categories)})")
    print(adata)
    print("=======================\n")

    return adata


if __name__=="__main__":
    from docopt import docopt

    args = docopt(__doc__)

    file = args["<file>"]
    name = args["--name"]
    labels = args["--labels"]
    out_file = args["--out-file"]

    adata = prep_dataset(name, file, labels)

    print(f"Writing data to '{out_file}'...")
    adata.write_h5ad(out_file)
