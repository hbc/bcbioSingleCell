#' Heatmap with Quantile Breaks
#'
#' @rdname quantileHeatmap
#' @name quantileHeatmap
#'
#' @details This is helpful for more usefully visualizing single cell data.
#'   Ideas and code from: http://slowkow.com/notes/heatmap-tutorial/
#'
#' @param object Matrix of data.
#' @param annotation Column annotations.
#' @param clusterRows Perform row clustering.
#' @param clusterCols Perform column clustering.
#'
#' @return `pheatmap()` object.
NULL



# Constructors ====
#' Create Breaks Based on Quantiles of the Data
#'
#' @param xs Numeric vector.
#' @param n The number of breaks to create.
#'
#' @return A vector of `n` quantile breaks.
.quantileBreaks <- function(xs, n = 10L) {
    breaks <- quantile(xs, probs = seq(0L, 1L, length.out = n))
    breaks[!duplicated(breaks)]
}



# Methods ====
#' @rdname quantileHeatmap
#' @export
setMethod("quantileHeatmap", "matrix", function(
    object,
    annotation = NA,
    clusterRows = TRUE,
    clusterCols = TRUE) {
    mat <- object
    if (isTRUE(quantile)) {
        matBreaks <- .quantileBreaks(mat)
    }
    if (isTRUE(clusterRows)) {
        matClusterRows <- dendsort(hclust(dist(mat)))
    }
    else {
        matClusterRows <- FALSE
    }
    if (isTRUE(clusterCols)) {
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
             color = inferno(length(matBreaks) - 1L),
             show_colnames = FALSE,
             show_rownames = FALSE)
})
