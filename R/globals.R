globalVariables(".")

packageVersion <- packageVersion("bcbioSingleCell")

#' @importFrom basejump.globals lanePattern
lanePattern <- basejump.globals::lanePattern

#' @importFrom bcbioBase projectDirPattern
projectDirPattern <- bcbioBase::projectDirPattern

# FIXME Take this out
separatorBar <- basejump.developer::separator()

requiredAssays <- "counts"

geom <- c("histogram", "ecdf", "violin", "ridgeline", "boxplot")
# c("violin", "ridgeline", "ecdf", "histogram", "boxplot")

# Trailing number is to match cellranger output.
barcodePattern <- ")_([ACGT_]{6,})(_[0-9]+)?$"

# `nCount` column is bcbioSingleCell class specific.
metricsCols <- c(
    "nUMI",
    "nGene",
    "nCoding",
    "nMito",
    "log10GenesPerUMI",
    "mitoRatio"
)
