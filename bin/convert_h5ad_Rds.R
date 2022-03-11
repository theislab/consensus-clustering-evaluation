#!/usr/bin/env Rscript

"
Convert H5AD to Rds

Usage:
    convert_h5ad_Rds.R --out-file=<path> [options] <file>

Options:
    -h --help             Show this screen.
    --out-file=<path>     Path to output file.
" -> doc

suppressPackageStartupMessages({
    library("zellkonverter")
})

if (sys.nframe() == 0) {
    args <- docopt::docopt(doc)

    file     <- args[["<file>"]]
    out_file <- args[["--out-file"]]

    reticulate::use_python(Sys.which("python"))
    anndata <- reticulate::import("anndata")
    message("Reading AnnData from '", file, "'...")
    adata <- anndata$read_h5ad(file)
    sce <- AnnData2SCE(adata, X_name = "counts", verbose = TRUE)
    message("Writing SingleCellExperiment to '", out_file, "'...")
    saveRDS(sce, out_file)
    message("Done!")
}
