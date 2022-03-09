import os.path
import urllib.request
import scanpy

url = "https://zenodo.org/record/5139561/files/COMBAT-CITESeq-DATA.h5ad?download=1"
full_file = "COMBAT_full.h5ad"
sample = "S00052-Ja005E-PBCa"
sample_col = "scRNASeq_sample_ID"
labels_col = "Annotation_cluster_name"
out_file = "COMBAT_S00052-Ja005E-PBCa.h5ad"

if not os.path.isfile(full_file):
    print(f"Downloading {full_file} from {url}...")
    urllib.request.urlretrieve(url, full_file)
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

print("Storing counts as X matrix...")
adata.X = adata.layers["raw"].copy()

print("Removing genes with 0 counts...")
scanpy.pp.filter_genes(adata, min_counts=1)

print("Final dataset:")
print(adata)

print(f"Saving dataset to {out_file}...")
adata.write_h5ad(out_file)

print("Done!")
