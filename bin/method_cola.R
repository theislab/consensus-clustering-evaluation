#!/usr/bin/env Rscript

"
Run the cola method

Usage:
    method_cola.R --out-file=<path> --labels=<str> [options] <file>

Options:
    -h --help             Show this screen.
    --out-file=<path>     Path to output file.
    --labels=<str>        Column of obs containing cell labels.
    --ncpus=<cpus>        Number of CPUs to use [default: 1].
" -> doc

suppressPackageStartupMessages({
    library("cola")
    library("SingleCellExperiment")
})

run_cola <- function(sce, labels, ncpus) {

    message("Normalising counts...")
    sce <- scuttle::logNormCounts(sce)

    message("Preparing matrix...")
    mat <- adjust_matrix(logcounts(sce))

    message("Calculating partitions with ", ncpus, " CPUS(s)...")
    n_labels <- length(unique(colData(sce)[[labels]]))
    min_k <- max(2, n_labels - 10)
    max_k <- n_labels + 10
    message("Min k: ", min_k, ", Max k: ", max_k)
    res <- consensus_partition(mat, k = seq(min_k, max_k), cores = ncpus)

    message("Selecting best k...")
    best_k <- suggest_best_k(res)
    message("Best k: ", best_k)

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
    ncpus     <- as.numeric(args[["--ncpus"]])

    message("Reading data from '", file, "'...")
    sce <- readRDS(file)
    sce <- run_cola(sce, labels, ncpus)
    message("Writing data to '", out_file, "'...")
    saveRDS(sce, out_file)
    message("Done!")
}
