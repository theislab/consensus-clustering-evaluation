import os
import urllib.request
import gzip
import shutil
import scanpy

url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE194122&format=file&file=GSE194122%5Fopenproblems%5Fneurips2021%5Fcite%5FBMMC%5Fprocessed%2Eh5ad%2Egz"
full_file = "Neurips_full.h5ad"
sample = "site2_donor1_cite"
sample_col = "Samplename"
labels_col = "cell_type"
out_file = "Neurips_site2-donor1-cite.h5ad"

if not os.path.isfile(full_file):
    print(f"Downloading {full_file} from {url}...")
    urllib.request.urlretrieve(url, full_file + ".gz")

    print("Decompressing downloaded file...")
    with gzip.open(full_file + ".gz", "rb") as f_compressed:
        with open(full_file, "wb") as f_decompressed:
            shutil.copyfileobj(f_compressed, f_decompressed)

    os.remove(full_file + ".gz")
else:
    print(f"{full_file} already exists")

print("Reading full dataset...")
adata = scanpy.read_h5ad(full_file)
print(adata)

print(f"Subsetting to {sample} sample...")
adata = adata[adata.obs[sample_col] == sample, :].copy()

print("Subsetting to labels with at least 20 cells...")
label_counts = adata.obs[labels_col].value_counts()
keep_labels = label_counts.index[label_counts.gt(19)]
adata = adata[adata.obs[labels_col].isin(keep_labels), :].copy()

print(f"Subsetting to GEX features...")
adata = adata[:, adata.var["feature_types"] == "GEX", ].copy()

print("Storing counts as X matrix...")
adata.X = adata.layers["counts"].copy()

print("Removing genes with 0 counts...")
scanpy.pp.filter_genes(adata, min_counts=1)

print("Final dataset:")
print(adata)

print(f"Saving dataset to {out_file}...")
adata.write_h5ad(out_file)

print("Done!")
