suppressPackageStartupMessages({
    library(splatter)      # v1.18.2
    library(zellkonverter) # v1.4.0
})

sim <- splatSimulateGroups(
    batchCells = 10000,
    group.prob = c(0.4, 0.4, 0.05, 0.05, 0.04, 0.03, 0.02, 0.01),
    nGenes     = 12000,
    seed       = 1
)

writeH5AD(
    sim,
    "sim-rare.h5ad",
    X_name   = "counts",
    assays   = FALSE,
    metadata = FALSE,
    colData  = c("Cell", "Group"),
    rowData  = "Gene",
    verbose  = TRUE
)
