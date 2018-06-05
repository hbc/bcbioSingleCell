globalVariables(".")

packageVersion <- packageVersion("bcbioSingleCell")

# Trailing number is to match cellranger output
barcodePattern <- ")_([ACGT_]{6,})(_[0-9]+)?$"

requiredAssays <- "counts"

# Empty sample metadata support (e.g. for splatter simulation SCE)
unknownSampleData <- data.frame(
    sampleID = "unknown",
    sampleName = "unknown",
    interestingGroups = "unknown",
    row.names = "unknown",
    stringsAsFactors = TRUE
)
