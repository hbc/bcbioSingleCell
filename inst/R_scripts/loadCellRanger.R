library(bcbioSingleCell)
bcb <- loadCellRanger(
    "data/cellranger",
    refDataDir = "annotations/refdata-cellranger-mm10-1.2.0",
    sampleMetadataFile = "meta/sample_metadata.xlsx",
    interestingGroups = c("genotype", "age"))
saveData(bcb)
