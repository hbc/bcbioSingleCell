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
plot_total_counts <- function(...) {
    .Deprecated("plot_read_counts")
    plot_read_counts(...)
}



#' @rdname deprecated
#' @export
plot_total_vs_detected <- function(...) {
    .Deprecated("plot_reads_vs_genes")
    plot_reads_vs_genes(...)
}
