.applyFilterCutoffs <- function(object) {
    .checkFilterParams(object)

    # Apply cell filtering cutoffs
    filterCells <- metadata(object)[["filterCells"]]
    if (is.null(filterCells)) {
        stop("'filterCells' metadata is 'NULL'", call. = FALSE)
    }
    cells <- intersect(
        colnames(object),
        metadata(object)[["filterCells"]]
    )
    if (is.null(cells)) {
        stop("'NULL' cells passed filtering", call. = FALSE)
    }
    object <- object[, cells]
    metadata(object)[["filterCells"]] <- cells

    # Apply gene filtering cutoffs
    # Warn here instead of stop, since we didn't define in earlier versions
    filterGenes <- metadata(object)[["filterGenes"]]
    if (is.null(filterGenes)) {
        warning("'filterGenes' metadata is 'NULL', call. = FALSE")
    }
    genes <- intersect(
        rownames(object),
        metadata(object)[["filterGenes"]]
    )
    if (is.null(genes)) {
        warning("'NULL' genes passed filtering", call. = FALSE)
    } else {
        object <- object[genes, ]
        metadata(object)[["filterGenes"]] <- genes
    }

    object
}
