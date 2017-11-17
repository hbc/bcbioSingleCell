.applyFilterCutoffs <- function(object) {
    .checkFilterParams(object)

    # Apply cell filtering cutoffs
    filterCells <- metadata(object)[["filterCells"]]
    if (is.null(filterCells)) {
        stop("'filterCells' metadata is 'NULL'", call. = FALSE)
    }
    cells <- intersect(
        colnames(object),
        metadata(object)[["filterCells"]]) %>%
        sort()
    if (is.null(cells)) {
        warning(paste(
            "'NULL' cells passed filtering.",
            "Resetting 'filterCells' metadata to all cells."
        ), call. = FALSE)
        cells <- colnames(object)
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
        metadata(object)[["filterGenes"]]) %>%
        sort()
    if (is.null(genes)) {
        warning(paste(
            "'NULL' genes passed filtering.",
            "Resetting 'filterGenes' metadata to include all genes."
        ), call. = FALSE)
        genes <- rownames(object)
    }
    object <- object[genes, ]
    metadata(object)[["filterGenes"]] <- genes

    object
}
