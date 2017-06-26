#' Deprecated functions
#'
#' @keywords internal
#' @author Michael Steinbaugh
#'
#' @param ... Passthrough parameters



#' @rdname deprecated
#' @export
barcode_metrics <- function(...) {
    .Deprecated("metrics")
    metrics(...)
}



#' @rdname deprecated
#' @export
plot_mito_counts <- function(...) {
    .Deprecated("plot_mito_ratio")
    plot_mito_ratio(...)
}



#' @rdname deprecated
#' @export
plot_total_cells <- function(...) {
    .Deprecated("plot_cell_counts")
    plot_cell_counts(...)
}



#' @rdname deprecated
#' @export
plot_total_vs_detected <- function(...) {
    .Deprecated("plot_umis_vs_genes")
    plot_umis_vs_genes(...)
}
