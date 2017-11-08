.checkFilterParams <- function(object) {
    filterParams <- metadata(object)[["filterParams"]]
    if (is.null(filterParams)) {
        stop(paste(
            "'filterCells()' must be run on the object"
        ), call. = FALSE)
    }
    if (!is.numeric(filterParams)) {
        stop("Filter parameters must all be numeric", call. = FALSE)
    }
}
