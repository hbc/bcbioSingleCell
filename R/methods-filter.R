# TODO Arrange plots into a grid?

#' Filter cellular barcodes
#'
#' Apply gene detection, mitochondrial abundance, and novelty score cutoffs to
#' cellular barcodes.
#'
#' @rdname filter
#' @author Michael Steinbaugh
#'
#' @param object [bcbioSCDataSet].
#' @param reads Minimum number of total read counts per cell.
#' @param min_genes Minimum number of genes detected.
#' @param max_genes Maximum number of genes detected.
#' @param mito_ratio Maximum relative mitochondrial abundance (`0-1` scale).
#' @param novelty Minimum novelty score.
#' @param show Show summary statistics and plots.
#'
#' @return New [bcbioSCDataSet].
#' @export
setMethod("filter", "bcbioSCDataSet", function(
    object,
    reads = 1000,
    min_genes = 500,
    max_genes = 5000,
    mito_ratio = 0.2,
    novelty = 0.8,
    show = TRUE) {
    name <- deparse(substitute(object))

    # Cellular barcode count
    counts(object) %>%
        ncol %>%
        paste("cellular barcodes detected") %>%
        message

    # Subset metrics
    metrics <- metrics(object) %>%
        as.data.frame %>%
        rownames_to_column %>%
        filter(.data[["total_counts"]] > !!reads,
               # FIXME Include option to filter by `coding_counts`?
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

    # Note that barcode identifiers are columns in the count matrix
    counts <- counts(object)[, rownames(metrics)]

    # Set up the filtered object
    # FIXME Have to reset the slots using replace methods
    filtered <- object
    filtered[["barcodes"]] <- NULL
    filtered[["metrics"]] <- metrics
    filtered[["counts"]] <- counts

    # Save filtering parameters
    filtered[["filter"]] <- list(
        reads = reads,
        min_genes = min_genes,
        max_genes = max_genes,
        mito_ratio = mito_ratio,
        novelty = novelty)

    if (isTRUE(show)) {
        writeLines(c(
            paste(name, "filtering parameters:"),
            paste0("- `> ", reads, "` total read counts per cell"),
            paste0("- `> ", min_genes, "` genes per cell"),
            paste0("- `< ", max_genes, "` genes per cell"),
            paste0("- `< ", mito_ratio, "` mitochondrial abundance ratio"),
            paste0("- `> ", novelty, "` novelty score")))
        plot_cell_counts(filtered)
        plot_read_counts(filtered, min = reads)
        plot_genes_detected(filtered, min = min_genes, max = max_genes)
        plot_reads_vs_genes(filtered) %>% show
        plot_mito_ratio(filtered, max = mito_ratio)
        plot_novelty(filtered, min = novelty)
    }

    filtered
})
