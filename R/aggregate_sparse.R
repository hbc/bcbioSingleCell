#' Aggregate sparse matrix operations
#'
#' @rdname aggregate_sparse
#' @author Rory Kirchner
#' @author Michael Steinbaugh
#'
#' @param sparse Sparse counts matrix.
#'
#' @return Modified sparse matrix.



#' @rdname aggregate_sparse
#' @description Collapse features with the same feature id by summing them.
#' @param featureids New feature IDs for the collapsed features.
#' @export
aggregate_sparse_features <- function(sparse, featureids) {
    rownames(sparse) <- featureids
    sparse <- sparse[!is.na(rownames(sparse)), ]
    aggregate.Matrix(sparse, rownames(sparse), fun = "sum")
}



#' @rdname aggregate_sparse
#' @description Pool cellular technical replicate counts.
#' @param cellids New cellular ids of the collapsed cells, one for each cell.
#' @export
aggregate_sparse_replicates <- function(sparse, cellids) {
    tsparse <- t(sparse)
    rownames(tsparse) <- cellids
    aggregate.Matrix(tsparse, cellids, fun = "sum") %>% t
}
