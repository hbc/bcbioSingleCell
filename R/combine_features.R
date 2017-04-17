#' Combine features of a sparse matrix
#'
#' This collapsed features with the same feature id by summing them
#'
#' @author Rory Kirchner
#'
#' @param sparse A sparse matrix of counts
#' @param featureids New feature ids for the collapsed features
#'
#' @return Sparse matrix
#' @export
combine_features <- function(sparse, featureids) {
  rownames(sparse) <- featureids
  sparse <- sparse[!is.na(rownames(sparse)), ]
  combine <- aggregate.Matrix(sparse, rownames(sparse), fun = "sum")
  return(combine)
}
