globalVariables(".")

packageVersion <- packageVersion("bcbioSingleCell")

lanePattern <- basejump::lanePattern
projectDirPattern <- bcbioBase::projectDirPattern
separatorBar <- basejump::separatorBar
updateMessage <- basejump::updateMessage

# Trailing number is to match cellranger output.
barcodePattern <- ")_([ACGT_]{6,})(_[0-9]+)?$"

requiredAssays <- "counts"

# `nCount` column is bcbioSingleCell class specific.
metricsCols <- c(
    "nUMI",
    "nGene",
    "nCoding",
    "nMito",
    "log10GenesPerUMI",
    "mitoRatio"
)

# Empty sample metadata support (e.g. for splatter simulation SCE).
unknownSampleData <- data.frame(
    sampleID = "unknown",
    sampleName = "unknown",
    interestingGroups = "unknown",
    row.names = "unknown",
    stringsAsFactors = TRUE
)
