suppressPackageStartupMessages({
    library(splatter)      # v1.18.2
    library(zellkonverter) # v1.4.0
})

mini_sim <- splatSimulateGroups(
    batchCells = 2000,
    nGenes     = 5000,
    group.prob = c(0.05, 0.3, 0.25, 0.4),
    seed       = 1
)

writeH5AD(
    mini_sim,
    "mini-sim.h5ad",
    X_name = "counts",
    assays = FALSE,
    metadata = FALSE,
    colData = c("Cell", "Group"),
    rowData = "Gene",
    verbose = TRUE
)
