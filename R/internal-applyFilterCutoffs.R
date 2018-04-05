.applyFilterCutoffs <- function(object) {
    validObject(object)

    params <- metadata(object)[["filterParams"]]
    if (!is.list(params)) {
        abort("`filterCells()` hasn't been applied to this dataset")
    }

    # Ensure all params are numeric
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
        x = colnames(object),
        y = metadata(object)[["filterCells"]]
    ) %>%
        sort()
    assert_is_non_empty(cells)
    object <- object[, cells]
    metadata(object)[["filterCells"]] <- cells

    # Apply gene filtering cutoffs =============================================
    filterGenes <- metadata(object)[["filterGenes"]]
    assert_is_character(filterGenes)
    assert_is_non_empty(filterGenes)
    genes <- intersect(
        x = rownames(object),
        y = metadata(object)[["filterGenes"]]
    ) %>%
        sort()
    assert_is_non_empty(genes, severity = "warning")
    if (is.null(genes)) {
        genes <- rownames(object)
    }
    object <- object[genes, ]
    metadata(object)[["filterGenes"]] <- genes

    object
}
