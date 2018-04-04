# nolint start
#
# Seurat markers per sample
# Michael Steinbaugh
# 2018-03-12
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

library(rmarkdown)
library(bcbioSingleCell)

data_dir <- file.path("data", Sys.Date())

# Sample subsets was saved in clustering_per_sample.R
loadData(subsets, dir = data_dir)

# Render R Markdown reports per bcbioSingleCell subset file
invisible(mapply(
    file = subsets,
    name = names(subsets),
    FUN = function(file, name) {
        seurat_name <- paste(name, "seurat", sep = "_")
        seurat_file <- file.path(data_dir, paste0(seurat_name, ".rda"))
        stopifnot(file.exists(seurat_file))
        render(
            input = "markers.Rmd",
            output_file = paste(seurat_name, "markers.html", sep = "_"),
            output_format = "html_document",
            params = list(seurat_file = seurat_file)
        )
    }
))
