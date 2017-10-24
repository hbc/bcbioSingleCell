library(pbapply)
library(rmarkdown)
source("setup.R")

# Enforce caching
opts_chunk[["set"]](
    audodep = TRUE,
    cache = TRUE,
    cache.lazy = FALSE)

loadData(bcb)

sampleIDs <- sampleMetadata(bcb) %>%
    pull("sampleID")
sampleSubsets <- pblapply(seq_along(sampleIDs), function(a) {
    sampleID <- sampleIDs[[a]]
    subset <- selectSamples(bcb, sampleID = sampleID)
    assignAndSaveData(name = sampleID, object = subset)
    sampleID
}) %>%
    unlist()
saveData(sampleSubsets)

# Render RMarkdown reports per bcbioSingleCell subset file
pblapply(seq_along(sampleSubsets), function(a) {
    bcbName <- sampleSubsets[[a]]
    bcbFile <- paste0("data/", bcbName, ".rda")
    seuratName <- paste(bcbName, "seurat", sep = "_")
    render(input = "seurat.Rmd",
           output_file = paste0(bcbName, ".html"),
           output_format = "html_document",
           params = list(
               bcbFile = bcbFile,
               seuratName = seuratName
           ))
})
