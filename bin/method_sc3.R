#!/usr/bin/env Rscript

"
Run the SC3 method

Usage:
    method_sc3.R --out-file=<path> [options] <file>

Options:
    -h --help             Show this screen.
    --out-file=<path>     Path to output file.
" -> doc

suppressPackageStartupMessages({
    library("SC3")
    library("SingleCellExperiment")
})

run_sc3 <- function(sce) {

    sce <- scuttle::logNormCounts(sce)
    rowData(sce)$feature_symbol <- rownames(sce)
    sce <- sc3_prepare(sce, n_cores = 1)
    sce <- sc3_estimate_k(sce)
    n_clusters <- metadata(sce)$sc3$k_estimation
    sce <- sc3_calc_dists(sce)
    sce <- sc3_calc_transfs(sce)
    sce <- sc3_kmeans(sce, ks = n_clusters)
    sce <- sc3_calc_consens(sce)

    colData(sce)$Cluster <- colData(sce)[paste0("sc3_", n_clusters, "_clusters")]

    return(sce)
}

if (sys.nframe() == 0) {
    args <- docopt::docopt(doc)

    file      <- args[["<file>"]]
    out_file  <- args[["--out-file"]]

    message("Reading data from '", file, "'...")
    sce <- readRDS(file)
    sce <- run_sc3(sce)
    message("Writing data to '", out_file, "'...")
    saveRDS(sce, out_file)
    message("Done!")
}
