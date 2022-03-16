#!/usr/bin/env Rscript

"
Run the SIMLR method

Usage:
    method_simlr.R --out-file=<path> [options] <file>

Options:
    -h --help             Show this screen.
    --out-file=<path>     Path to output file.
    --labels=<str>        Column of obs containing cell labels [default: Label].
    --ncpus=<cpus>        Number of CPUs to use [default: 1].
" -> doc

suppressPackageStartupMessages({
    library("SIMLR")
    library("SingleCellExperiment")
})

run_simlr <- function(sce, labels, ncpus) {

    message("Normalising counts...")
    sce <- scuttle::logNormCounts(sce)

    cores_ratio = ncpus / parallel::detectCores()
    message("Running SIMLR with ", ncpus, " CPU(s), cores ratio: ", cores_ratio)

    mat <- as.matrix(logcounts(sce))
    message("Estimating number of clusters...")
    if (ncol(sce) > 1000) {
        message("Using 1000 randomly selected cells to estimate k")
        set.seed(1)
        sel_mat <- mat[, sample(ncol(sce), 1000)]
    } else {
        sel_mat <- mat
    }
    n_labels <- length(unique(colData(sce)[[labels]]))
    min_k <- max(2, n_labels - 5)
    max_k <- n_labels + 5
    message("Min k: ", min_k, ", Max k: ", max_k)
    k_ests <- SIMLR_Estimate_Number_of_Clusters(
        sel_mat,
        NUMC        = seq(min_k, max_k),
        cores.ratio = cores_ratio
    )
    message("Selecting estimated k...")
    min_k1 <- which.min(k_ests$K1) + 1
    min_k2 <- which.min(k_ests$K2) + 1
    k_est <- max(min_k1, min_k2)
    message("Estimated k: ", k_est)

    message("Running SIMLR...")
    results <- SIMLR_Large_Scale(mat, k_est)

    data.frame(
        Cell    = colnames(sce),
        Label   = colData(sce)$Label,
        Cluster = results$y$cluster
    )
}

if (sys.nframe() == 0) {
    args <- docopt::docopt(doc)

    file      <- args[["<file>"]]
    out_file  <- args[["--out-file"]]
    labels    <- args[["--labels"]]
    ncpus     <- as.numeric(args[["--ncpus"]])

    message("Reading data from '", file, "'...")
    sce <- readRDS(file)
    print(sce)
    clusters <- run_simlr(sce, labels, ncpus)
    message("Writing clusters to '", out_file, "'...")
    write.table(
        clusters,
        file      = out_file,
        quote     = FALSE,
        sep       = "\t",
        row.names = FALSE
    )
    message("Done!")
}
