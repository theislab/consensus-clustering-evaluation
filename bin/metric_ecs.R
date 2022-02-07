#!/usr/bin/env Rscript

"
Run the Element-Centric Clustering Similarity (ECS) metric

Usage:
    metric_ecs.R --out-file=<path> --dataset=<str> --labels=<str> --method=<str> [options] <file>

Options:
    -h --help             Show this screen.
    --out-file=<path>     Path to output file.
    --dataset=<str>       Name of the dataset.
    --labels=<str>        Column of obs containing cell labels.
    --method=<str>        Name of the method.
" -> doc

suppressPackageStartupMessages({
    library("SingleCellExperiment")
})

run_ecs <- function(sce, dataset, labels, method) {

    message("Calculating ECS...")
    score <- ClustAssess::element_sim(sce[[labels]], sce[["Cluster"]])

    data.frame(
        Dataset = dataset,
        Method  = method,
        Metric  = "ECS",
        Score   = score
    )
}

if (sys.nframe() == 0) {
    args <- docopt::docopt(doc)

    file      <- args[["<file>"]]
    out_file  <- args[["--out-file"]]
    dataset   <- args[["--dataset"]]
    labels    <- args[["--labels"]]
    method    <- args[["--method"]]

    message("Reading data from '", file, "'...")
    sce <- readRDS(file)
    message("Read data:")
    print(sce)
    results <- run_ecs(sce, dataset, labels, method)
    message("Writing data to '", out_file, "'...")
    write.table(results, out_file, quote = FALSE, sep = "\t", row.names = FALSE)
    message("Done!")
}
