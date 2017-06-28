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
    # Assays ====
    # Subset the cellular barcodes in contained in the count matrices by
    # colData. When packaging a bcbioSCDataSet, very low stringency subsetting
    # of the raw data is applied, based on the cellular barcode summary metrics,
    # to save disk space. This is very helpful when loading unfiltered counts,
    # such as a raw dataset from 10X Cell Ranger. By default, the bcbio-nextgen
    # pipeline applies cellular barcode filtering server-side, so this step
    # should have little to no effect.
    sparse_counts <- sparse_counts[, rownames(col_data)]
    if (!identical(colnames(sparse_counts), rownames(col_data))) {
        stop("colData mismatch")
    }


    # rowData ====
    # Subset the annotable by the genes present in the counts matrix
    row_data <- row_data[rownames(sparse_counts), ] %>%
        set_rownames(rownames(sparse_counts))


    # Metadata ====
    if (is.null(metadata)) {
        metadata <- SimpleList()
    }

    # Update automatic metadata slots
    metadata[["date"]] <- Sys.Date()
    metadata[["wd"]] <- getwd()
    metadata[["hpc"]] <- detect_hpc()
    metadata[["session_info"]] <- sessionInfo()

    # Check for retired Ensembl identifiers, which can happen when a more recent
    # annotable build is used than the genome build. If present, store these
    # identifiers in the metadata.
    if (any(is.na(row_data[["ensgene"]]))) {
        warning("Ensembl identifier degradation detected")
        metadata[["retired_ensgene"]] <- row_data %>%
            as_tibble %>%
            filter(is.na(.data[["ensgene"]])) %>%
            pull("rowname") %>%
            sort
    }


    # Return ====
    if (!identical(colnames(sparse_counts), rownames(col_data))) {
        stop("colData mismatch")
    }
    if (!identical(rownames(sparse_counts), rownames(row_data))) {
        stop("rowData mismatch")
    }
    SummarizedExperiment(
        assays = SimpleList(
            sparse_counts = sparse_counts),
        colData = as(col_data, "DataFrame"),
        rowData = as(row_data, "DataFrame"),
        metadata = as(metadata, "SimpleList"))
}



.row_data <- function(object) {
    rowData(object) %>% set_rownames(rownames(object))
}
