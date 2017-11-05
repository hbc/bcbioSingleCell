.applyFilterCutoffs <- function(object) {
    .checkFilterCells(object)
    cells <- intersect(
        colnames(object),
        metadata(object)[["filterCells"]])
    if (!is.null(cells)) {
        object <- object[, cells]
    } else {
        warning("'NULL' cells passed filtering", call. = FALSE)
    }
    genes <- intersect(
        rownames(object),
        metadata(object)[["filterGenes"]])
    if (!is.null(genes)) {
        object <- object[genes, ]
    } else {
        warning("'NULL' genes passed filtering", call. = FALSE)
    }
    object
}
