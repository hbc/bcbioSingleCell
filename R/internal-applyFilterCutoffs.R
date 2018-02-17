.applyFilterCutoffs <- function(object) {
    .checkFilterParams(object)

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
