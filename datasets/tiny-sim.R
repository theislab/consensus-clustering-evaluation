suppressPackageStartupMessages({
    library(splatter)      # v1.18.2
    library(zellkonverter) # v1.4.0
})

sim <- splatSimulateGroups(
    batchCells = 500,
    nGenes     = 800,
    group.prob = c(0.5, 0.5),
    seed       = 1,
    lib.loc    = 9.5
)

writeH5AD(
    sim,
    "tiny-sim.h5ad",
    X_name   = "counts",
    assays   = FALSE,
    metadata = FALSE,
    colData  = c("Cell", "Group"),
    rowData  = "Gene",
    verbose  = TRUE
)
