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

    file      <- args[["<file>"]]
    out_file  <- args[["--out-file"]]

    sce <- readH5AD(file, verbose = TRUE)
    message("Writing data to '", out_file, "'...")
    saveRDS(sce, out_file)
    message("Done!")
}
