cacheURL <- "http://bcbiosinglecell.seq.cloud"
mapply(
    FUN = function(cacheURL, file) {
        # Download file to testthat directory
        if (!file_exists(file)) {
            download.file(url = paste(cacheURL, file, sep = "/"), destfile = file)
        }
        # Load R Data file
        if (grepl("\\.rda", file)) {
            load(file)
        }
    },
    file = c(
        "bcb.rda",
        "filtered.rda",
        "known_markers_detected.rda",
        "pooled.rda",
        "seurat.rda",
        "seurat_all_markers.rda",
        "seurat_all_markers_original.rda",
        "top.rda"
    ),
    MoreArgs = list(cacheURL = cacheURL)
)
