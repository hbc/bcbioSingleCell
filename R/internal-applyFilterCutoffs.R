.applyFilterCutoffs <- function(object) {
    .checkFilterParams(object)

    # Apply cell filtering cutoffs =============================================
    filterCells <- metadata(object)[["filterCells"]]
    if (is.null(filterCells)) {
        abort("`filterCells` metadata is `NULL`")
    }

    cells <- intersect(
        colnames(object),
        metadata(object)[["filterCells"]]) %>%
        sort()
    if (!length(cells)) {
        warn(paste(
            "No cells passed filtering.",
            "Resetting `filterCells` metadata to all cells."
        ))
        cells <- colnames(object)
    }

    object <- object[, cells]
    metadata(object)[["filterCells"]] <- cells

    # Apply gene filtering cutoffs =============================================
    # Warn here instead of stop, since we didn't define in earlier versions
    filterGenes <- metadata(object)[["filterGenes"]]
    if (is.null(filterGenes)) {
        warn("`filterGenes` metadata is `NULL`")
    }

    genes <- intersect(
        rownames(object),
        metadata(object)[["filterGenes"]]) %>%
        sort()
    if (is.null(genes)) {
        warn(paste(
            "`NULL` genes passed filtering.",
            "Resetting `filterGenes` metadata to include all genes."
        ))
        genes <- rownames(object)
    }

    object <- object[genes, ]
    metadata(object)[["filterGenes"]] <- genes

    object
}
