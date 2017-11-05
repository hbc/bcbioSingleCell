.checkFilterCells <- function(object) {
    if (is.null(metadata(object)[["filterParams"]])) {
        stop(paste(
            "Run 'filterCells()' on the bcbioSingleCell dataset",
            "before selecting samples"
        ), call. = FALSE)
    }
    if (is.null(metadata(object)[["filterCells"]])) {
        stop("'filterCells' metadata is 'NULL'", call. = FALSE)
    }
    # Warn here instead of stop, since we didn't set in earlier versions
    if (is.null(metadata(object)[["filterGenes"]])) {
        warning("'filterGenes' metadata is 'NULL', call. = FALSE")
    }
}
