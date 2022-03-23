suppressPackageStartupMessages({
    library(splatter)      # v1.18.2
    library(zellkonverter) # v1.4.0
})

sim <- splatSimulate(
    batchCells = 10000,
    nGenes     = 12000,
    seed       = 1
)

sim$Group <- "Group1"

writeH5AD(
    sim,
    "sim-blob.h5ad",
    X_name   = "counts",
    assays   = FALSE,
    metadata = FALSE,
    colData  = c("Cell", "Group"),
    rowData  = "Gene",
    verbose  = TRUE
)
