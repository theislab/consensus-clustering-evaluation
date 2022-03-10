#!/usr/bin/env Rscript

"
Convert Rds to H5AD

Usage:
    convert_Rds_h5ad.R --out-file=<path> [options] <file>

Options:
    -h --help             Show this screen.
    --out-file=<path>     Path to output file.
" -> doc

suppressPackageStartupMessages({
    library("zellkonverter")
})

if (sys.nframe() == 0) {
    args <- docopt::docopt(doc)

    file      <- args[["<file>"]]
    out_file  <- args[["--out-file"]]

    reticulate::use_python(Sys.which("python"))
    anndata <- reticulate::import("anndata")
    message("Reading SingleCellExperiment from '", file, "'...")
    sce <- readRDS(file)
    adata <- SCE2AnnData(sce, verbose = TRUE)
    message("Writing AnnData to '", out_file, "'...")
    adata$write_h5ad(out_file)
    message("Done!")
}
