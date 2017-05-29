combine_features <- function(...) {
    .Deprecated("aggregate_sparse_features")
    aggregate_sparse_features(...)
}

combine_technical_replicates <- function(...) {
    .Deprecated("aggregate_sparse_replicates")
    aggregate_sparse_replicates(...)
}

filter_barcodes <- function(...) {
    .Deprecated("filter_cellular_barcodes")
    filter_cellular_barcodes(...)
}

load_csv_with_rownames <- function(...) {
    .Deprecated("read_csv_with_rownames")
    read_csv_with_rownames(...)
}

load_sparsecounts <- function(...) {
    .Deprecated("read_counts")
    read_counts(...)
}

make_bcbio_skeleton <- function(...) {
    .Deprecated("create_bcbio_skeleton")
    create_bcbio_skeleton(...)
}

plot_mito_counts <- function(...) {
    .Deprecated("plot_mito")
    plot_mito(...)
}

plot_mito_counts_boxplot <- function(...) {
    .Deprecated("plot_mito_boxplot")
    plot_mito_boxplot(...)
}

plot_mito_counts_histogram <- function(...) {
    .Deprecated("plot_mito_histogram")
    plot_mito_histogram(...)
}

plot_mito_counts_scatterplot <- function(...) {
    .Deprecated("plot_mito_scatterplot")
    plot_mito_scatterplot(...)
}

plot_total_cells <- function(...) {
    .Deprecated("plot_cell_counts")
    plot_cell_counts(...)
}

plot_total_vs_detected <- function(...) {
    .Deprecated("plot_reads_vs_genes")
    plot_reads_vs_genes(...)
}
