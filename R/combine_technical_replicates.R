#' Combine cellular technical replicates in a sparse matrix
#'
#' @author Rory Kirchner
#'
#' @param sparse a sparse matrix of counts
#' @param cellids new cellular ids of the collapsed cells, one for each cell
#'
#' @return sparse matrix with technical replicates combined
#' @export
combine_technical_replicates <- function(sparse, cellids) {
    tsparse <- t(sparse)
    rownames(tsparse) <- cellids
    return(t(aggregate.Matrix(tsparse, cellids, fun = "sum")))
}
