#' Aggregation methods
#'
#' @rdname aggregate
#'
#' @author Rory Kirchner, Michael Steinbaugh
#'
#' @param object [bcbioSCDataSet].
#'
#' @return [bcbioSCDataSet].
#' @export
setMethod("aggregate_replicates", "bcbioSCDataSet", function(object) {
    stop("Draft function to be added in a future update")
})



#' Collapse features with the same feature id by summing them.
#'
#' @rdname aggregate
#'
#' @param sparse Sparse counts matrix.
#' @param featureids New feature IDs for the collapsed features.
#'
#' @return Aggregated sparse matrix.
.aggregate_features <- function(sparse, featureids) {
    rownames(sparse) <- featureids
    sparse <- sparse[!is.na(rownames(sparse)), ]
    aggregate.Matrix(sparse, rownames(sparse), fun = "sum")
}



#' Pool cellular technical replicate counts.
#'
#' @rdname aggregate
#'
#' @param cellids New cellular ids of the collapsed cells, one for each cell.
#'
#' @return Aggregated sparse matrix.
.aggregate_replicates <- function(sparse, cellids) {
    tsparse <- t(sparse)
    rownames(tsparse) <- cellids
    aggregate.Matrix(tsparse, cellids, fun = "sum") %>% t
}
