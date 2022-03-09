#!/usr/bin/env Rscript

"
Run the SIMLR method

Usage:
    method_simlr.R --out-file=<path> [options] <file>

Options:
    -h --help             Show this screen.
    --out-file=<path>     Path to output file.
    --ncpus=<cpus>        Number of CPUs to use [default: 1].
" -> doc

suppressPackageStartupMessages({
    library("SIMLR")
    library("SingleCellExperiment")
})

run_simlr <- function(sce, ncpus) {

    message("Normalising counts...")
    sce <- scuttle::logNormCounts(sce)

    cores_ratio = ncpus / parallel::detectCores()
    message("Running SIMLR with ", ncpus, " CPU(s), cores ratio: ", cores_ratio)

    message("Estimating number of clusters...")
    k_ests <- SIMLR_Estimate_Number_of_Clusters(
        logcounts(sce),
        NUMC        = 2:30,
        cores.ratio = cores_ratio
    )
    min_k1 <- which.min(k_ests$K1) + 1
    min_k2 <- which.min(k_ests$K2) + 1
    k_est <- max(min_k1, min_k2)

    message("Running SIMLR...")
    results <- SIMLR_Large_Scale(logcounts(sce), k_est)
    colData(sce)$Cluster <- results$y$cluster

    return(sce)
}

if (sys.nframe() == 0) {
    args <- docopt::docopt(doc)

    file      <- args[["<file>"]]
    out_file  <- args[["--out-file"]]
    ncpus     <- as.numeric(args[["--ncpus"]])

    message("Reading data from '", file, "'...")
    sce <- readRDS(file)
    sce <- run_simlr(sce, ncpus)
    message("Writing data to '", out_file, "'...")
    saveRDS(sce, out_file)
    message("Done!")
}
