library(bcbioSingleCell)
library(pbapply)
library(rmarkdown)

prepareSingleCellTemplate()
source("setup.R")

loadData(bcb, dir = "data")

sampleIDs <- sampleMetadata(bcb) %>%
    pull(sampleID)
sampleSubsets <- pblapply(seq_along(sampleIDs), function(a) {
    sampleID <- sampleIDs[[a]]
    subset <- selectSamples(bcb, sampleID = sampleID)
    assignAndSaveData(name = sampleID, object = subset, dir = "data")
    sampleID
}) %>%
    unlist()
saveData(sampleSubsets)

# Render RMarkdown reports per bcbioSingleCell subset file
pblapply(seq_along(sampleSubsets), function(a) {
    bcbName <- sampleSubsets[[a]]
    bcbFile <- file.path("data", paste0(bcbName, ".rda"))
    seuratName <- paste(bcbName, "seurat", sep = "_")
    render(input = "clustering.Rmd",
           output_file = paste0(bcbName, ".html"),
           output_format = "html_document",
           params = list(
               bcbFile = bcbFile,
               seuratName = seuratName
           ))
}) %>%
    invisible()
