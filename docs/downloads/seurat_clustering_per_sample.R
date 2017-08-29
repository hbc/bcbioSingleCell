library(bcbioSingleCell)

# Generate per sample SCSubset objects
data(bcbFiltered)
meta <- sampleMetadata(bcbFiltered)

sampleIDs <- pull(meta, "sampleID")
sampleNames <- pull(meta, "sampleName")
camelSampleNames <- camel(sampleNames)

sampleSubsets <- pbsapply(seq_along(sampleIDs), function(a) {
    sampleID <- sampleIDs[[a]]
    sampleName <- sampleNames[[a]]
    camelSampleName <- camelSampleNames[[a]]
    subset <- selectSamples(
        bcbFiltered,
        sampleID = sampleID)
    # Check for enough cells
    if (dim(subset)[[2L]] < 200L) {
        warning(paste(
            sampleName,
            "sample has < 200 cells...",
            "skipping"),
            call. = FALSE)
        return(NA)
    }
    assignAndSaveData(camelSampleName, subset, env = globalenv())
    c(sampleID = camelSampleName)
}) %>%
    na.omit %>%
    as.character
rm(bcbFiltered)
saveData(sampleSubsets)

# Render RMarkdown reports per SCSubset
message("Rendering RMarkdown reports")
pbsapply(seq_along(sampleSubsets), function(a) {
    subsetName <- sampleSubsets[[a]]
    message(subsetName)
    render(input = "seurat.Rmd",
           output_file = paste0(subsetName, ".html"),
           output_format = "html_document",
           params = list(subsetName = subsetName))
})
