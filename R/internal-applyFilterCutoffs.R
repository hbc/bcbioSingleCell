.applyFilterCutoffs <- function(object) {
    validObject(object)
    params <- metadata(object)[["filterParams"]]

    # Warn and early return
    if (is.null(params)) {
        warn("Use `filterCells()` apply filtering cutoffs")
        return(object)
    }

    # TODO Require user to `updateObject()` for this in the future
    # `filterParams` metadata was stored as a named numeric vector up until
    # v0.0.28. We changed to storing as a list in v0.0.29, to allow for per
    # sample cutoffs.
    if (is.numeric(params)) {
        params <- as.list(params)
    }
    assert_is_list(params)

    # Ensure all params are numeric
    # TODO Switch to assertive method here
    if (!all(vapply(
        X = params,
        FUN = is.numeric,
        FUN.VALUE = logical(1L)
    ))) {
        abort("Filter parameters must be numeric")
    }

    # Apply cell filtering cutoffs =============================================
    filterCells <- metadata(object)[["filterCells"]]
    assert_is_character(filterCells)
    assert_is_non_empty(filterCells)
    cells <- intersect(
        colnames(object),
        metadata(object)[["filterCells"]]) %>%
        sort()
    assert_is_non_empty(cells)
    object <- object[, cells]
    metadata(object)[["filterCells"]] <- cells

    # Apply gene filtering cutoffs =============================================
    filterGenes <- metadata(object)[["filterGenes"]]
    # Warn here instead of stop, since we didn't define in earlier versions
    assert_is_character(filterCells, severity = "warning")
    assert_is_non_empty(filterGenes, severity = "warning")
    genes <- intersect(
        rownames(object),
        metadata(object)[["filterGenes"]]) %>%
        sort()
    assert_is_non_empty(genes, severity = "warning")
    if (is.null(genes)) {
        genes <- rownames(object)
    }
    object <- object[genes, ]
    metadata(object)[["filterGenes"]] <- genes

    object
}
