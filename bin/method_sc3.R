#!/usr/bin/env Rscript

"
Run the SC3 method

Usage:
    method_sc3.R --out-file=<path> [options] <file>

Options:
    -h --help             Show this screen.
    --out-file=<path>     Path to output file.
    --ncpus=<cpus>        Number of CPUs to use [default: 1]
" -> doc

suppressPackageStartupMessages({
    library("SC3")
    library("SingleCellExperiment")
})

run_sc3 <- function(sce, ncpus) {

    message("Calculating log normalised counts...")
    sce <- scuttle::logNormCounts(sce)
    rowData(sce)$feature_symbol <- rownames(sce)
    counts(sce) <- as.matrix(counts(sce))
    logcounts(sce) <- as.matrix(logcounts(sce))
    message("Running SC3 with ", ncpus, " CPU(s)")
    sce <- sc3_prepare(sce, n_cores = ncpus, rand_seed = 1)
    sce <- sc3_estimate_k(sce)
    n_clusters <- metadata(sce)$sc3$k_estimation
    sce <- sc3_calc_dists(sce)
    sce <- sc3_calc_transfs(sce)
    sce <- sc3_kmeans(sce, ks = n_clusters)
    sce <- sc3_calc_consens(sce)

    if (nrow(sce) > 5000) {
        message("Transferring clusters using SVM...")
        sce <- sc3_run_svm(sce, ks = n_clusters)
    }

    data.frame(
        Cell    = colnames(sce),
        Label   = colData(sce)$Label,
        Cluster = colData(sce)[[paste0("sc3_", n_clusters, "_clusters")]]
    )
}

if (sys.nframe() == 0) {
    args <- docopt::docopt(doc)

    file      <- args[["<file>"]]
    out_file  <- args[["--out-file"]]
    ncpus     <- as.numeric(args[["--ncpus"]])

    message("Reading data from '", file, "'...")
    sce <- readRDS(file)
    print(sce)
    clusters <- run_sc3(sce, ncpus)
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
