# nolint start
#
# Seurat markers per sample
# Michael Steinbaugh
# 2017-12-14
#
# Clustering must be performed before running this script!
#
# Latest version of this script is available here:
# script <- system.file(
#     file.path("R_scripts", "markers_per_sample.R"),
#     package = "bcbioSingleCell")
# file.edit(script)
#
# nolint end

library(bcbioSingleCell)
library(pbapply)
library(rmarkdown)

dataDir <- "data"

# sampleSubsets was saved in clustering_per_sample.R
loadData(sampleSubsets, dir = dataDir)

# Render R Markdown reports per bcbioSingleCell subset file
pblapply(seq_along(sampleSubsets), function(a) {
    bcbName <- sampleSubsets[[a]]
    seuratName <- paste(bcbName, "seurat", sep = "_")
    seuratFile <- file.path(dataDir, paste0(seuratName, ".rda"))
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
