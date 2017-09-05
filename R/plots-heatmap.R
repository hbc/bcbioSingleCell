
##' heatmap with quantile breaks
##'
##' This is helpful for more usefully visualizing single cell data.
##' ideas and code from: http://slowkow.com/notes/heatmap-tutorial/
##' @title heatmap with quantile breaks
##'
##' @param mat matrix of data
##' @param annotation column annotation
##' @param clusterRows perform row clustering
##' @param clusterCols perform column clustering
##' @importFrom pheatmap pheatmap
##' @importFrom viridis inferno
##' @importFrom dendsort dendsort
##' @return pheatmap object
##' @export
##' @author Rory Kirchner
quantileHeatmap <- function(mat, annotation=NA, clusterRows=TRUE,
                        clusterCols=TRUE) {
  if(isTRUE(quantile)) {
    matBreaks = .quantileBreaks(mat)
  }
  if(isTRUE(cluster_rows)) {
    matClusterRows <- dendsort(hclust(dist(mat)))
  }
  else {
    matClusterRows <- FALSE
  }
  if(isTRUE(cluster_cols)) {
    matClusterCols <- dendsort(hclust(dist(t(mat))))
  }
  else {
    matClusterCols <- FALSE
  }
  pheatmap(mat,
           annotation_col = annotation,
           cluster_cols = matClusterCols,
           cluster_rows = matClusterRows,
           breaks = matBreaks,
           color = inferno(length(matBreaks) - 1),
           show_colnames = FALSE,
           show_rownames = FALSE)
}

##' create breaks based on quantiles of the data
##'
##' @title quantile breaks
##' @param xs a vector of numbers
##' @param n the number of breaks to create
##' @return a vector of n quantile breaks
##' @author Rory Kirchner
.quantileBreaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}
