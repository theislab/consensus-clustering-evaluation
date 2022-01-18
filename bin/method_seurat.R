#!/usr/bin/env Rscript

"
Run the Seurat method

Usage:
    method_seurat.R --out-file=<path> [options] <file>

Options:
    -h --help             Show this screen.
    --out-file=<path>     Path to output file.
" -> doc

suppressPackageStartupMessages({
    library("Seurat")
    library("SingleCellExperiment")
})

run_seurat <- function(seurat) {

    message("Scaling and log transforming...")
    seurat <- NormalizeData(seurat, scale.factor = 10000)
    message("Finding highly variable genes...")
    seurat <- FindVariableFeatures(seurat)
    message("Found ", length(VariableFeatures(seurat)), " HVGs")
    message("Scaling data...")
    seurat <- ScaleData(seurat)
    message("Calculating PCA...")
    seurat <- RunPCA(seurat)
    message("Building SNN graph...")
    seurat <- FindNeighbors(seurat)
    message("Running Louvain community detection...")
    seurat <- FindClusters(seurat)

}

if (sys.nframe() == 0) {
    args <- docopt::docopt(doc)

    file      <- args[["<file>"]]
    out_file  <- args[["--out-file"]]

    message("Reading data from '", file, "'...")
    sce <- readRDS(file)
    message("Converting to Seurat '", file, "'...")
    seurat <- as.Seurat(sce, data = NULL)
    seurat <- run_seurat(seurat)
    message("Seurat to SingleCellExperiment '", file, "'...")
    sce <- as.SingleCellExperiment(seurat)
    colData(sce)$Cluster <- colData(sce)$ident
    message("Writing data to '", out_file, "'...")
    saveRDS(sce, out_file)
    message("Done!")
}
