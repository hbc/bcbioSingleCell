cacheURL <- "http://bcbiosinglecell.seq.cloud"
files <- c(
    "bcb.rda",
    "filtered.rda",
    "known_markers_detected.rda",
    "pooled.rda",
    "seurat.rda",
    "seurat_all_markers.rda",
    "seurat_all_markers_original.rda",
    "top.rda")
mapply(
    FUN = function(cacheURL, file, envir) {
        # Download file to testthat directory
        if (!file.exists(file)) {
            download.file(
                url = paste(cacheURL, file, sep = "/"),
                destfile = file)
        }
        # Load R Data file
        if (grepl("\\.rda$", file)) {
            message(paste("Loading", file))
            load(file, envir = envir)
        }
    },
    file = files,
    MoreArgs = list(cacheURL = cacheURL, envir = environment())
)
