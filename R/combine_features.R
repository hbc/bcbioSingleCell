##' combine features of a sparse matrix
##'
##' this collapsed features with the same feature id by summing them
##'
##' @param sparse a sparse matrix of counts
##' @param featureids new feature ids for the collapsed features
##' @return sparse matrix
##' @author Rory Kirchner
##' @export
combine_features <- function(sparse, featureids) {
  rownames(sparse) <- featureids
  sparse <- sparse[!is.na(rownames(sparse)), ]
  return(aggregate.Matrix(sparse, featureids, fun="sum"))
}
