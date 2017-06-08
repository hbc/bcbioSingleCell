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
combine_features <- function(...) {
    .Deprecated("aggregate_sparse_features")
    aggregate_sparse_features(...)
}



#' @rdname deprecated
#' @export
combine_technical_replicates <- function(...) {
    .Deprecated("aggregate_sparse_replicates")
    aggregate_sparse_replicates(...)
}



#' @rdname deprecated
#' @export
filter_barcodes <- function(...) {
    .Deprecated("filter_cellular_barcodes")
    filter_cellular_barcodes(...)
}



#' @rdname deprecated
#' @export
load_csv_with_rownames <- function(...) {
    .Deprecated("read_csv_with_rownames")
    read_csv_with_rownames(...)
}



#' @rdname deprecated
#' @export
load_sparsecounts <- function(...) {
    .Deprecated("read_counts")
    read_counts(...)
}



#' @rdname deprecated
#' @export
make_bcbio_skeleton <- function(...) {
    .Deprecated("create_bcbio_skeleton")
    create_bcbio_skeleton(...)
}



#' @rdname deprecated
#' @export
plot_mito_counts <- function(...) {
    .Deprecated("plot_mito_ratio")
    plot_mito_ratio(...)
}



#' @rdname deprecated
#' @export
plot_mito_counts_boxplot <- function(...) {
    .Deprecated("plot_mito_ratio_boxplot")
    plot_mito_ratio_boxplot(...)
}



#' @rdname deprecated
#' @export
plot_mito_counts_histogram <- function(...) {
    .Deprecated("plot_mito_ratio_histogram")
    plot_mito_ratio_histogram(...)
}



#' @rdname deprecated
#' @export
plot_mito_counts_scatterplot <- function(...) {
    .Deprecated("plot_mito_ratio_scatterplot")
    plot_mito_ratio_scatterplot(...)
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



#' @rdname deprecated
#' @export
read_10x_counts <- function(...) {
    .Deprecated("read_10x")
    read_10x(...)
}
