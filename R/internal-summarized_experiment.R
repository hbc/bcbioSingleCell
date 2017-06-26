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


    # colData ====
    # Pre-filter the sparse counts matrix, based on barcode metrics (colData).
    # Very low stringency, requires detailed follow-up quality control analysis.
    # At this point, we're just removing zero count barcodes.
    sparse_counts <- sparse_counts[, rownames(col_data)]
    identical(colnames(sparse_counts), rownames(col_data))


    # rowData ====
    row_data <- row_data[rownames(sparse_counts), ] %>%
        set_rownames(rownames(sparse_counts))
    # Check for retired Ensembl identifiers, which can happen when a more recent
    # annotable build is used than the genome build. This can happen when
    # importing 10X Cell Ranger counts.
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
