#' Aggregation methods
#'
#' @rdname aggregate_replicates
#' @keywords internal
#'
#' @author Rory Kirchner, Michael Steinbaugh
#'
#' @param object Primary object.
#'
#' @return [bcbioSCDataSet].
#' @export
setMethod("aggregate_replicates", "bcbioSCDataSet", function(object) {
    stop("Draft function to be added in a future update")
})



#' Collapse features with the same feature id by summing them.
#'
#' @rdname aggregate_replicates
#' @usage NULL
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
#' @rdname aggregate_replicates
#' @usage NULL
#'
#' @param cellids New cellular ids of the collapsed cells, one for each cell.
#'
#' @return Aggregated sparse matrix.
.aggregate_replicates <- function(sparse, cellids) {
    tsparse <- t(sparse)
    rownames(tsparse) <- cellids
    aggregate.Matrix(tsparse, cellids, fun = "sum") %>% t
}
