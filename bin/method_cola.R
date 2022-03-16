#!/usr/bin/env Rscript

"
Run the cola method

Usage:
    method_cola.R --out-file=<path> [options] <file>

Options:
    -h --help             Show this screen.
    --out-file=<path>     Path to output file.
    --labels=<str>        Column of obs containing cell labels [default: Label].
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
    mat <- adjust_matrix(as.matrix(logcounts(sce)))

    message("Running cola with ", ncpus, " CPUS(s)...")
    if (ncol(mat) > 1000) {
        message("Using 1000 randomly selected cells to estimate k")
        set.seed(1)
        sel_mat <- mat[, sample(ncol(mat), 1000)]
    } else {
        sel_mat <- mat
    }
    n_labels <- length(unique(colData(sce)[[labels]]))
    min_k <- max(2, n_labels - 5)
    max_k <- n_labels + 5
    message("Min k: ", min_k, ", Max k: ", max_k)
    res <- consensus_partition(sel_mat, k = seq(min_k, max_k), cores = ncpus)

    message("Selecting best k...")
    best_k <- suggest_best_k(res)
    message("Best k: ", best_k)

    message("Performing final partitioning...")
    res <- consensus_partition(mat, k = best_k, cores = ncpus)

    message("Assigning clusters...")
    classes <- get_classes(res)
    colData(sce)$Cluster <- classes$class

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
