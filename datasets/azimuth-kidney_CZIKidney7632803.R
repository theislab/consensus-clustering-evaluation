suppressPackageStartupMessages({
    library(Seurat)               # v4.1.0
    library(SingleCellExperiment) # v1.16.0
    library(zellkonverter)        # v1.4.0
})

url        <- "https://seurat.nygenome.org/azimuth/demo_datasets/kidney_demo_stewart.rds"
full_file  <- "azimuth-kidney_full.Rds"
sample     <- "CZIKidney7632803"
sample_col <- "orig.ident"
labels_col <- "celltype"
out_file   <- "azimuth-kidney_CZIKidney7632803.h5ad"

if (!fs::file_exists(full_file)) {
    message(glue::glue("Downloading {full_file} from {url}..."))
    download.file(url, full_file)
} else {
    message(glue::glue("{full_file} already exists"))
}

message("Reading full dataset...")
seurat <- readRDS(full_file)
print(seurat)

message("Converting to SingleCellExperiment...")
sce <- as.SingleCellExperiment(seurat)

message(glue::glue("Subsetting to {sample} sample..."))
sce <- sce[, colData(sce)[[sample_col]] == sample]

print("Subsetting to labels with at least 20 cells...")
label_counts <- table(colData(sce)[[labels_col]])
keep_labels <- names(label_counts)[label_counts >= 20]
sce <- sce[, colData(sce)[[labels_col]] %in% keep_labels]

print("Removing genes with 0 counts...")
total_counts <- rowSums(counts(sce))
sce <- sce[total_counts > 0, ]

writeH5AD(
    sce,
    out_file,
    X_name   = "counts",
    assays   = FALSE,
    metadata = FALSE,
    colData  = c(sample_col, labels_col),
    verbose  = TRUE
)
