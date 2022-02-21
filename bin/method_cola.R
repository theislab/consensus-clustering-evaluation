#!/usr/bin/env Rscript

"
Run the cola method

Usage:
    method_cola.R --out-file=<path> --labels=<str> [options] <file>

Options:
    -h --help             Show this screen.
    --out-file=<path>     Path to output file.
    --labels=<str>        Column of obs containing cell labels.
" -> doc

suppressPackageStartupMessages({
    library("cola")
    library("SingleCellExperiment")
})

run_cola <- function(sce, labels) {

    message("Normalising counts...")
    sce <- scuttle::logNormCounts(sce)

    message("Preparing matrix...")
    mat <- adjust_matrix(logcounts(sce))

    message("Calculating partitions...")
    n_labels <- length(unique(colData(sce)[[labels]]))
    min_k <- max(2, n_labels - 5)
    max_k <- n_labels + 5
    res <- consensus_partition(mat, k = seq(min_k, max_k))

    message("Selecing best k...")
    best_k <- suggest_best_k(res)

    message("Assigning clusters...")
    classes <- get_classes(res)
    colData(sce)$Cluster <- classes[[paste0("k=", best_k)]]

    return(sce)
}

if (sys.nframe() == 0) {
    args <- docopt::docopt(doc)

    file      <- args[["<file>"]]
    out_file  <- args[["--out-file"]]
    labels    <- args[["--labels"]]

    message("Reading data from '", file, "'...")
    sce <- readRDS(file)
    sce <- run_cola(sce, labels)
    message("Writing data to '", out_file, "'...")
    saveRDS(sce, out_file)
    message("Done!")
}
