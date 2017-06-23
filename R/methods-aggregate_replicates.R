#' Sparse matrix aggregation operations
#'
#' @rdname sparse
#' @author Rory Kirchner, Michael Steinbaugh
#'
#' @param sparse Sparse counts matrix.
#'
#' @return ["dgCMatrix-class"].
#'
#' @description Collapse features with the same feature id by summing them.
#' @param featureids New feature IDs for the collapsed features.
.aggregate_features <- function(sparse, featureids) {
    rownames(sparse) <- featureids
    sparse <- sparse[!is.na(rownames(sparse)), ]
    aggregate.Matrix(sparse, rownames(sparse), fun = "sum")
}



#' @rdname sparse
#' @description Pool cellular technical replicate counts.
#' @param cellids New cellular ids of the collapsed cells, one for each cell.
.aggregate_replicates <- function(sparse, cellids) {
    tsparse <- t(sparse)
    rownames(tsparse) <- cellids
    aggregate.Matrix(tsparse, cellids, fun = "sum") %>% t
}



#' Aggregate replicates
#'
#' @rdname aggregate_replicates
#' @export
#'
#' @return New [bcbioSCDataSet].
setMethod("aggregate_replicates", "bcbioSCDataSet", function(object) {
    stop("Draft function")
})
