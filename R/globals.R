globalVariables(".")

packageVersion <- packageVersion("bcbioSingleCell")

lanePattern <- basejump::lanePattern
projectDirPattern <- bcbioBase::projectDirPattern
separatorBar <- basejump::separatorBar
updateMessage <- basejump::updateMessage

requiredAssays <- "counts"

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

# We're grep matching against these camel case variants here to automatically
# sanitize `colData()` into sample-level `sampleData()`.
clusterCols <- c(
    "^ident$",
    "^origIdent$",
    "res[.0-9]^",
    "^sScore$",
    "^g2mScore$",
    "^phase$"
)

# Empty sample metadata support (e.g. for splatter simulation SCE).
unknownSampleData <- data.frame(
    sampleID = "unknown",
    sampleName = "unknown",
    interestingGroups = "unknown",
    row.names = "unknown",
    stringsAsFactors = TRUE
)
