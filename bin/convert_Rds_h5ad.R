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

    message("Writing file from '", file, "'...")
    sce <- readRDS(file)
    writeH5AD(sce, out_file, verbose = TRUE)
    message("Done!")
}
