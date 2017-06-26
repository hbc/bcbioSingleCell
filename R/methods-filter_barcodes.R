# TODO Arrange plots into a grid?

#' Filter cellular barcodes
#'
#' Apply gene detection, mitochondrial abundance, and novelty score cutoffs to
#' cellular barcodes.
#'
#' @rdname filter_barcodes
#' @author Michael Steinbaugh
#'
#' @param object [bcbioSCDataSet].
#' @param ... Additional parameters.
#'
#' @param reads Minimum number of total read counts per cell.
#' @param min_genes Minimum number of genes detected.
#' @param max_genes Maximum number of genes detected.
#' @param mito_ratio Maximum relative mitochondrial abundance (`0-1` scale).
#' @param novelty Minimum novelty score.
#' @param show Show summary statistics and plots.
#'
#' @return Filtered [bcbioSCDataSet], with low quality cellular barcodes that
#'   don't pass quality control cutoffs removed.
#' @export
setMethod("filter_barcodes", "bcbioSCDataSet", function(
    object,
    reads = 1000,
    min_genes = 500,
    max_genes = 10000,
    mito_ratio = 0.1,
    novelty = 0.8,
    show = TRUE) {
    name <- deparse(substitute(object))
    sparse_counts <- counts(object)

    # Cellular barcode count
    ncol(sparse_counts) %>%
        paste("cellular barcodes in dataset") %>%
        message

    # Subset metrics
    metrics <- metrics(object) %>%
        rownames_to_column %>%
        filter(.data[["total_counts"]] > !!reads,
               .data[["genes_detected"]] > !!min_genes,
               .data[["genes_detected"]] < !!max_genes,
               .data[["mito_ratio"]] < !!mito_ratio,
               .data[["log10_detected_per_count"]] > !!novelty) %>%
        column_to_rownames

    # Subset counts
    if (!nrow(metrics)) {
        stop("Barcode filtering was too stringent. No cells remaining.")
    }
    message(paste(nrow(metrics), "cellular barcodes passed filtering"))

    # Sparse
    sparse_size_original <- object.size(sparse_counts) %>%
        format(units = "auto")
    sparse_counts <- sparse_counts[, rownames(metrics)]
    sparse_size_filtered <- object.size(sparse_counts) %>%
        format(units = "auto")
    message(paste("Sparse matrix size shrunk from",
                  sparse_size_original, "to", sparse_size_filtered))

    # colData ====
    col_data <- colData(object) %>%
        .[rownames(metrics), ] %>%
        DataFrame

    # rowData ====
    row_data <- rowData(object) %>%
        as.data.frame %>%
        set_rownames(names(object)) %>%
        DataFrame

    # Metadata ====
    meta <- metadata(object)
    meta[["filtering_criteria"]] <- SimpleList(
        reads = reads,
        min_genes = min_genes,
        max_genes = max_genes,
        mito_ratio = mito_ratio,
        novelty = novelty)

    # SummarizedExperiment ====
    se <- .summarized_experiment(
        sparse_counts = sparse_counts,
        col_data = col_data,
        row_data = row_data,
        metadata = meta)

    # bcbioSCDataSet ====
    bcb <- new("bcbioSCDataSet", se)

    message(paste(
        "bcbioSCDataSet:",
        format(object.size(bcb), units = "auto")
    ))

    # Show summary statistics and plots, if desired
    if (isTRUE(show)) {
        writeLines(c(
            paste(name, "filtering parameters:"),
            paste0("- `> ", reads, "` total read counts per cell"),
            paste0("- `> ", min_genes, "` genes per cell"),
            paste0("- `< ", max_genes, "` genes per cell"),
            paste0("- `< ", mito_ratio, "` mitochondrial abundance ratio"),
            paste0("- `> ", novelty, "` novelty score")))
        plot_cell_counts(bcb)
        plot_read_counts(bcb, min = reads)
        plot_genes_detected(bcb, min = min_genes, max = max_genes)
        plot_reads_vs_genes(bcb)
        plot_mito_ratio(bcb, max = mito_ratio)
        plot_novelty(bcb, min = novelty)
    }

    bcb
})
