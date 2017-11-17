.checkFilterParams <- function(object) {
    filterParams <- metadata(object)[["filterParams"]]
    if (is.null(filterParams)) {
        stop(paste(
            "'filterCells()' hasn't been applied to this dataset"
        ), call. = FALSE)
    }
    if (!is.numeric(filterParams)) {
        stop("Filter parameters must all be numeric", call. = FALSE)
    }
}
