#' Deprecated Functions
#'
#' @rdname deprecated
#' @name deprecated
#' @keywords internal
#'
#' @return Varies by function.
NULL



# 0.0.14 ====
#' @rdname deprecated
#' @export
barcode_metrics <- function() {
    .Deprecated("metrics")
}

#' @rdname deprecated
#' @export
filter_barcodes <- function() {
    .Deprecated("filterCells")
}

#' @rdname deprecated
#' @export
plot_genes_detected <- function() {
    .Deprecated("plotGenesPerCell")
}

#' @rdname deprecated
#' @export
plot_mito_counts <- function() {
    .Deprecated("plotMitoRatio")
}

#' @rdname deprecated
#' @export
plot_total_cells <- function() {
    .Deprecated("plotCellCounts")
}

#' @rdname deprecated
#' @export
plot_total_vs_detected <- function() {
    .Deprecated("plotUMIsVsGenes")
}
