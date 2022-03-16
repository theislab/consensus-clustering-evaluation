suppressPackageStartupMessages({
    library(splatter)      # v1.18.2
    library(zellkonverter) # v1.4.0
})

sim <- splatSimulatePaths(
    batchCells = 10000,
    group.prob = c(0.5, 0.5),
    nGenes     = 12000,
    seed       = 1
)

writeH5AD(
    sim,
    "sim-path.h5ad",
    X_name   = "counts",
    assays   = FALSE,
    metadata = FALSE,
    colData  = c("Cell", "Group"),
    rowData  = "Gene",
    verbose  = TRUE
)
