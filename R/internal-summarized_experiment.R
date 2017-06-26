#' Package scRNA-Seq counts into [SummarizedExperiment]
#'
#' @rdname summarized_experiment
#' @keywords internal
#'
#' @author Michael Steinbaugh
#'
#' @param sparse_counts Gene-level sparse counts matrix.
#' @param col_data Cellular barcode metrics.
#' @param row_data Annotable.
#' @param metadata Custom metadata.
#'
#' @return [SummarizedExperiment].
.summarized_experiment <- function(
    sparse_counts,
    col_data,
    row_data,
    metadata = NULL) {
    message("Packaging SummarizedExperiment")
    # Coerce metadata to [SimpleList], if necessary
    if (!is.null(metadata) & class(metadata) != "SimpleList") {
        metadata <- as(metadata, "SimpleList")
    }

    # Pre-filter the sparse counts matrix, based on barcode metrics (colData).
    # Very low stringency, requires detailed follow-up quality control analysis.
    # At this point, we're just removing cellular barcodes with very low
    # read counts or gene detection.
    sparse_counts <- sparse_counts[, rownames(col_data)]
    if (!identical(colnames(sparse_counts), rownames(col_data))) {
        stop("colData mismatch")
    }


    # rowData ====
    # Subset the annotable by the genes present in the counts matrix
    row_data <- row_data[rownames(sparse_counts), ] %>%
        set_rownames(rownames(sparse_counts))
    if (!identical(rownames(sparse_counts), rownames(row_data))) {
        stop("rowData mismatch")
    }

    # Check for retired Ensembl identifiers, which can happen when a more recent
    # annotable build is used than the genome build. If present, store these
    # identifiers in the metadata.
    if (any(is.na(row_data[["ensgene"]]))) {
        warning("Ensembl identifier degradation detected")
        metadata[["retired_ensgene"]] <- row_data %>%
            as.data.frame %>%
            rownames_to_column %>%
            filter(is.na(.data[["ensgene"]])) %>%
            pull("rowname") %>%
            sort
    }


    # Return [SummarizedExperiment] ====
    SummarizedExperiment(
        assays = SimpleList(
            sparse_counts = sparse_counts),
        colData = col_data,
        rowData = row_data,
        metadata = metadata)
}
