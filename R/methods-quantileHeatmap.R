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
#' @return `pheatmap()`.
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



.quantileHeatmap <- function(
    object,
    annotation = NA,
    clusterRows = TRUE,
    clusterCols = TRUE) {
    mat <- as.matrix(object)
    matBreaks <- .quantileBreaks(mat)
    # Dendrogram sorting can take a long time on large datasets
    if (isTRUE(clusterRows)) {
        matClusterRows <- dendsort(hclust(dist(mat)))
    } else {
        matClusterRows <- FALSE
    }
    if (isTRUE(clusterCols)) {
        matClusterCols <- dendsort(hclust(dist(t(mat))))
    } else {
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
}



# Methods ====
#' @rdname quantileHeatmap
#' @export
setMethod("quantileHeatmap", "dgCMatrix", .quantileHeatmap)



#' @rdname quantileHeatmap
#' @export
setMethod("quantileHeatmap", "matrix", .quantileHeatmap)



#' @rdname quantileHeatmap
#' @export
setMethod("quantileHeatmap", "seurat", function(
    object,
    annotation = NA,
    clusterRows = TRUE,
    clusterCols = TRUE) {
    # Use the raw counts
    .quantileHeatmap(
        object@raw.data,
        annotation = annotation,
        clusterRows = clusterRows,
        clusterCols = clusterCols)
})
