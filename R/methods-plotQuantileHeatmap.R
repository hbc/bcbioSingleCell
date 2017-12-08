#' Plot Heatmap with Quantile Breaks
#'
#' @details This is helpful for more usefully visualizing single cell data.
#'   Ideas and code from: http://slowkow.com/notes/heatmap-tutorial/
#'
#' @rdname plotQuantileHeatmap
#' @name plotQuantileHeatmap
#' @author Rory Kirchner
#'
#' @inheritParams AllGenerics
#'
#' @param object Matrix of data.
#' @param n The number of breaks to create.
#' @param annotationCol Column annotations [data.frame].
#' @param clusterRows Perform row clustering.
#' @param clusterCols Perform column clustering.
#' @param color Color palette function.
#'
#' @return [pheatmap::pheatmap()].
#'
#' @examples
#' load(system.file(
#'     file.path("inst", "extdata", "filtered.rda"),
#'     package = "bcbioSingleCell"))
#'
#' # Reduce the number of counts to speed up the working example
#' subset <- filtered[1:250, 1:250]
#'
#' # bcbioSingleCell
#' plotQuantileHeatmap(
#'     subset,
#'     n = 20,
#'     clusterRows = FALSE,
#'     clusterCols = FALSE)
NULL



# Constructors ====
#' Create Breaks Based on Quantiles of the Data
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom stats quantile
#'
#' @param x Numeric vector.
#' @param n The number of breaks to create.
#' @param unique Only return unique quantiles.
#'
#' @return A vector of `n` quantile breaks.
.quantileBreaks <- function(x, n = 10, unique = TRUE) {
    q <- quantile(x, probs = seq(0, 1, length.out = n))
    if (isTRUE(unique)) {
        q <- q[!duplicated(q)]
    }
    q
}



#' Quantile Heatmap Constructor
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom dendsort dendsort
#' @importFrom pheatmap pheatmap
#' @importFrom stats dist hclust
#' @importFrom viridis viridis
.plotQuantileHeatmap <- function(
    object,
    n = 10,
    annotationCol = NA,
    clusterRows = FALSE,
    clusterCols = FALSE,
    color = viridis::viridis) {
    if (!is.function(color)) {
        stop("'color' argument must contain a color palette function",
             call. = FALSE)
    }
    mat <- as.matrix(object)
    breaks <- .quantileBreaks(mat, n = n)

    # Dendrogram sorting can take a long time on large datasets
    if (isTRUE(clusterRows)) {
        clusterRows <- dendsort(hclust(dist(mat)))
    } else {
        clusterRows <- FALSE
    }
    if (isTRUE(clusterCols)) {
        clusterCols <- dendsort(hclust(dist(t(mat))))
    } else {
        clusterCols <- FALSE
    }

    pheatmap(
        mat,
        annotation_col = annotationCol,
        cluster_cols = clusterCols,
        cluster_rows = clusterRows,
        breaks = breaks,
        color = color(length(breaks)),
        show_colnames = FALSE,
        show_rownames = FALSE)
}



# Methods ====
#' @rdname plotQuantileHeatmap
#' @export
setMethod(
    "plotQuantileHeatmap",
    signature("bcbioSingleCell"),
    function(
        object,
        n = 10,
        annotationCol = NA,
        clusterRows = FALSE,
        clusterCols = FALSE) {
        counts <- counts(object)
        .plotQuantileHeatmap(
            object = counts,
            n = n,
            annotationCol = annotationCol,
            clusterRows = clusterRows,
            clusterCols = clusterCols)
    })



#' @rdname plotQuantileHeatmap
#' @export
setMethod(
    "plotQuantileHeatmap",
    signature("dgCMatrix"),
    .plotQuantileHeatmap)



#' @rdname plotQuantileHeatmap
#' @export
setMethod(
    "plotQuantileHeatmap",
    signature("matrix"),
    .plotQuantileHeatmap)



#' @rdname plotQuantileHeatmap
#' @export
setMethod(
    "plotQuantileHeatmap",
    signature("seurat"),
    function(
        object,
        n = 10,
        annotationCol = NA,
        clusterRows = FALSE,
        clusterCols = FALSE) {
        # Use the raw counts
        counts <- counts(object, normalized = FALSE)
        .plotQuantileHeatmap(
            object = counts,
            n = n,
            annotationCol = annotationCol,
            clusterRows = clusterRows,
            clusterCols = clusterCols)
    })
