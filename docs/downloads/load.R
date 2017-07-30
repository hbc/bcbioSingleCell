library(bcbioSinglecell)
bcb <- loadSinglecell(
    uploadDir = file.path("data", "final"),
    sampleMetadataFile = file.path("meta", "sample_metadata.xlsx"),
    interestingGroups = c("genotype", "treatment"),
    experimentName = "",
    researcher = "",
    principalInvestigator = "",
    author = getOption("author"),
    email = getOption("email"))
saveData(bcb)
