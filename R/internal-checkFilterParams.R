.checkFilterParams <- function(object) {
    filterParams <- metadata(object)[["filterParams"]]
    if (is.null(filterParams)) {
        abort("`filterCells()` hasn't been applied to this dataset")
    }
    if (!is.numeric(filterParams)) {
        abort("Filter parameters must all be numeric")
    }
}
