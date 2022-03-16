suppressPackageStartupMessages({
    library(splatter)      # v1.18.2
    library(zellkonverter) # v1.4.0
})

sim <- splatSimulateGroups(
    batchCells = 10000,
    group.prob = c(0.3, 0.2, 0.2, 0.1, 0.1, 0.1),
    nGenes     = 12000,
    seed       = 1
)

writeH5AD(
    sim,
    "sim-groups.h5ad",
    X_name   = "counts",
    assays   = FALSE,
    metadata = FALSE,
    colData  = c("Cell", "Group"),
    rowData  = "Gene",
    verbose  = TRUE
)
