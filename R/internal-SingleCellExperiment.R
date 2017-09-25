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
    metadata) {
    # Assays ====
    assays <- as.list(assays)
    assay <- assays[[1]]

    # Row data ====
    rowData <- as.data.frame(rowData) %>%
        .[rownames(assay), , drop = FALSE]
    # Check for gene mismatch
    if (!all(rownames(assay) %in% rownames(rowData))) {
        missing <- setdiff(rownames(assay), rownames(rowData))
        warning(paste(
            "rowData mismatch with assay slot:",
            paste0(toString(missing), "."),
            "These identifiers are missing in the current Ensembl release."
        ))
    }

    # Column data ====
    colData <- as.data.frame(colData) %>%
        .[colnames(assay), , drop = FALSE]

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
