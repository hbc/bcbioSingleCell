# Latest version of this script is available here:
# script <- system.file(
#     file.path("R_scripts", "markers_per_sample.R"),
#     package = "bcbioSingleCell")
# file.edit(script)

# Clustering must be performed before running this script!

library(bcbioSingleCell)
library(pbapply)
library(rmarkdown)

prepareSingleCellTemplate()
source("setup.R")
dataDir <- "data"

# sampleSubsets was saved in clustering_per_sample.R
loadData(sampleSubsets, dir = "data")

# Render RMarkdown reports per bcbioSingleCell subset file
pblapply(seq_along(sampleSubsets), function(a) {
    bcbName <- sampleSubsets[[a]]
    seuratName <- paste(bcbName, "seurat", sep = "_")
    seuratFile <- file.path("data", paste0(seuratName, ".rda"))
    if (!file.exists(seuratFile)) {
        stop(paste(
            "Seurat file missing:", seuratFile
        ), call. = FALSE)
    }
    render(input = "markers.Rmd",
           output_file = paste0(seuratName, "_markers.html"),
           output_format = "html_document",
           params = list(
               seuratFile = seuratFile
           ))
}) %>%
    invisible()
