library(bcbioSingleCell)
bcb <- loadSingleCell(
    uploadDir = "bcbio_indrop_rnaseq/final",
    interestingGroups = c("genotype", "treatment"),
    sampleMetadataFile = "meta/sample_metadata.xlsx",
    organism = "Homo sapiens",
    ensemblVersion = 90L
)
saveData(bcb, dir = file.path("data", Sys.Date()))
