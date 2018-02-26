cacheURL <- "http://bcbiosinglecell.seq.cloud"
mapply(
    FUN = function(cacheURL, file) {
        if (!file_exists(file)) {
            download.file(url = paste(cacheURL, file, sep = "/"), destfile = file)
        }
    },
    file = c(
        "bcb.rda",
        "filtered.rda",
        "knownMarkersDetected.rda",
        "pooled.rda",
        "seurat.rda",
        "seuratAllMarkers.rda",
        "seuratAllMarkersOriginal.rda",
        "topMarkers.rda"
    ),
    MoreArgs = list(cacheURL = cacheURL)
)
