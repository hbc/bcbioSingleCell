# nolint start
#
# Seurat Clustering per Sample
#
# Michael Steinbaugh
# 2018-02-26
#
# Latest version of this script is available here:
# script <- system.file(
#     "R_scripts/clustering_per_sample.R",
#     package = "bcbioSingleCell")
# file.edit(script)
#
# nolint end

library(rmarkdown)
library(bcbioSingleCell)

data_dir <- path("data", Sys.Date())

loadDataAsName(bcb = "bcb_filtered", dir = data_dir)
metadata(bcb)[["filterParams"]]

subsets <- subsetPerSample(bcb, dir = data_dir)
saveData(subsets, dir = data_dir)

# Render R Markdown reports per bcbioSingleCell subset file
invisible(mapply(
    file = subsets,
    name = names(subsets),
    FUN = function(file, name) {
        seurat_name <- paste(name, "seurat", sep = "_")
        render(input = "clustering.Rmd",
               output_file = paste0(name, ".html"),
               output_format = "html_document",
               params = list(
                   bcb_file = file,
                   seurat_name = seurat_name
               ))
    }
))
