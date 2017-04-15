
##' combine technical replicates
##'
##' @param counts a sparse matrix of counts
##' @param cellids new cellular ids of the collapsed cells
##' @return sparse matrix with technical replicates combined
##' @author Rory Kirchner
##' @export
combine_technical_replicates = function(counts, cellids) {
  tcounts = t(counts)
  rownames(tcounts) = cellids
  return(t(aggregate.Matrix(tcounts, cellids, fun="sum")))
}
