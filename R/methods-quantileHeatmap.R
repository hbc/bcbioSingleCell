#' Heatmap with Quantile Breaks
#'
#' @details This is helpful for more usefully visualizing single cell data.
#'   Ideas and code from: http://slowkow.com/notes/heatmap-tutorial/
#'
#' @rdname quantileHeatmap
#' @name quantileHeatmap
#' @author Rory Kirchner
#'
#' @inheritParams AllGenerics
#'
#' @param object Matrix of data.
#' @param annotation Column annotations.
#' @param clusterRows Perform row clustering.
#' @param clusterCols Perform column clustering.
#'
#' @return [pheatmap::pheatmap()].
NULL



# Constructors ====
#' Create Breaks Based on Quantiles of the Data
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom stats quantile
#'
#' @param xs Numeric vector.
#' @param n The number of breaks to create.
#'
#' @return A vector of `n` quantile breaks.
.quantileBreaks <- function(xs, n = 10) {
    xs %>%
        quantile(probs = seq(0, 1, length.out = n)) %>%
        .[!duplicated(.)]
}



#' Quantile Heatmap Constructor
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom dendsort dendsort
#' @importFrom pheatmap pheatmap
#' @importFrom stats dist hclust
#' @importFrom viridis inferno
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
    pheatmap(
        mat,
        annotation_col = annotation,
        cluster_cols = matClusterCols,
        cluster_rows = matClusterRows,
        breaks = matBreaks,
        color = inferno(length(matBreaks) - 1),
        show_colnames = FALSE,
        show_rownames = FALSE)
}



# Methods ====
#' @rdname quantileHeatmap
#' @export
setMethod(
    "quantileHeatmap",
    signature("dgCMatrix"),
    quantileHeatmap)



#' @rdname quantileHeatmap
#' @export
setMethod(
    "quantileHeatmap",
    signature("matrix"),
    .quantileHeatmap)



#' @rdname quantileHeatmap
#' @export
setMethod(
    "quantileHeatmap",
    signature("seurat"),
    function(
        object,
        annotation = NA,
        clusterRows = TRUE,
        clusterCols = TRUE) {
        # Use the raw counts
        .quantileHeatmap(
            slot(object, "raw.data"),
            annotation = annotation,
            clusterRows = clusterRows,
            clusterCols = clusterCols)
    })
