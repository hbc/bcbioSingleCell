# nolint start
#
# Seurat clustering per sample
# Michael Steinbaugh
# 2017-12-01
#
# Latest version of this script is available here:
# script <- system.file(
#     file.path("R_scripts", "clustering_per_sample.R"),
#     package = "bcbioSingleCell")
# file.edit(script)
#
# nolint end

library(bcbioSingleCell)
library(pbapply)
library(rmarkdown)

prepareSingleCellTemplate()
source("setup.R")
dataDir <- "data"

loadDataAsName(bcb = "bcb_filtered", dir = dataDir)
metadata(bcb)$filterParams

sampleNames <- sampleMetadata(bcb) %>%
    pull(sampleName)
sampleSubsets <- pblapply(seq_along(sampleNames), function(a) {
    sampleName <- as.character(sampleNames[[a]])
    sampleFileName <- snake(sampleName)
    subset <- selectSamples(bcb, sampleName = sampleName)
    assignAndSaveData(name = sampleFileName, object = subset, dir = dataDir)
    sampleFileName
}) %>%
    unlist()
saveData(sampleSubsets, dir = dataDir)

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
