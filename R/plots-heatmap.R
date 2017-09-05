
##' heatmap with quantile breaks
##'
##' This is helpful for more usefully visualizing single cell data.
##' ideas and code from: http://slowkow.com/notes/heatmap-tutorial/
##' @title heatmap with quantile breaks
##'
##' @param mat matrix of data
##' @param annotation column annotation
##' @param cluster_rows perform row clustering
##' @param cluster_cols perform column clustering
##' @importFrom pheatmap pheatmap
##' @importFrom viridis inferno
##' @importFrom dendsort dendsort
##' @return pheatmap object
##' @export
##' @author Rory Kirchner
quantileHeatmap <- function(mat, annotation=NA, cluster_rows=TRUE,
                        cluster_cols=TRUE) {
  if(isTRUE(quantile)) {
    mat_breaks = .quantile_breaks(mat)
  }
  if(isTRUE(cluster_rows)) {
    mat_cluster_rows = dendsort(hclust(dist(mat)))
  }
  else {
    mat_cluster_rows = FALSE
  }
  if(isTRUE(cluster_cols)) {
    mat_cluster_cols = dendsort(hclust(dist(t(mat))))
  }
  else {
    mat_cluster_cols = FALSE
  }
  pheatmap(mat,
           annotation_col=annotation,
           cluster_cols=mat_cluster_cols,
           cluster_rows=mat_cluster_rows,
           breaks=mat_breaks,
           color=inferno(length(mat_breaks) - 1),
           show_colnames=FALSE,
           show_rownames=FALSE)
}

##' create breaks based on quantiles of the data
##'
##' @title quantile breaks
##' @param xs a vector of numbers
##' @param n the number of breaks to create
##' @return a vector of n quantile breaks
##' @author Rory Kirchner
.quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}
