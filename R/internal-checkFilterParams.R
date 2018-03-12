.checkFilterParams <- function(object) {
    params <- metadata(object)[["filterParams"]]

    # Warn and early return on unfiltered dataset
    if (is.null(params)) {
        return(warn("`filterCells()` hasn't been applied to this dataset"))
    }

    # `filterParams` metadata was stored as a named numeric vector up until
    # v0.0.28. We changed to storing as a list in v0.0.29, to allow for per
    # sample cutoffs.
    if (is.numeric(params)) {
        params <- as.list(params)
    }

    # Ensure all params are numeric
    if (!all(vapply(
        X = params,
        FUN = is.numeric,
        FUN.VALUE = logical(1L)
    ))) {
        abort("Filter parameters must be numeric")
    }
}
