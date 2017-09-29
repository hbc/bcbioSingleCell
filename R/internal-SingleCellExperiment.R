#' Prepare SingleCellExperiment Object
#'
#' @inherit basejump::prepareSummarizedExperiment
#'
#' @seealso [basejump::prepareSummarizedExperiment()].
#'
#' @return [SingleCellExperiment].
#' @noRd
.SingleCellExperiment <- function(  # nolint
    assays,
    rowData,
    colData,
    metadata,
    subset = FALSE) {
    # Assays ====
    assays <- as.list(assays)
    assay <- assays[[1]]

    # Row data ====
    # Warn on any gene mismatches
    if (!all(rownames(assay) %in% rownames(rowData))) {
        missing <- setdiff(rownames(assay), rownames(rowData))
        warning(paste(
            "'rowData' mismatch with 'assay':",
            toString(missing)
        ), call. = FALSE)
    }
    rowData <- rowData %>%
        as.data.frame() %>%
        .[rownames(assay), , drop = FALSE] %>%
        set_rownames(rownames(assay))

    # Column data ====
    # Stop on any cell mismatches
    if (!all(colnames(assay) %in% rownames(colData))) {
        missing <- setdiff(colnames(assay), rownames(colData))
        stop(paste(
            "'colData' mismatch with 'assay':",
            toString(missing)
        ), call. = FALSE)
    }
    colData <- colData %>%
        as.data.frame() %>%
        .[colnames(assay), , drop = FALSE] %>%
        set_rownames(colnames(assay))

    # Metadata ====
    metadata <- as.list(metadata)
    # R session information
    metadata[["date"]] <- Sys.Date()
    metadata[["wd"]] <- getwd()
    metadata[["devtoolsSessionInfo"]] <-
        devtools::session_info(include_base = TRUE)
    metadata[["utilsSessionInfo"]] <-
        utils::sessionInfo()

    # Return ====
    SingleCellExperiment(
        assays = assays,
        rowData = rowData,
        colData = colData,
        metadata = metadata)
}
