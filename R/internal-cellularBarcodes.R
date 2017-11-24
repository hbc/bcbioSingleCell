#' Bind Cellular Barcodes
#'
#' @author Michael Steinbaugh
#' @keywords internal
#' @noRd
#'
#' @importFrom dplyr bind_rows group_by mutate
#' @importFrom parallel mclapply
#' @importFrom rlang !!
#' @importFrom tibble column_to_rownames
#'
#' @param list List of cellular barcodes.
#'
#' @return [data.frame].
.bindCellularBarcodes <- function(list) {
    list <- mclapply(seq_along(list), function(a) {
        # Add the sampleID as a column
        sampleID <- names(list)[[a]]
        list[[a]] %>%
            mutate(sampleID = !!sampleID)
    })
    df <- list %>%
        bind_rows() %>%
        as.data.frame() %>%
        mutate(
            rowname = paste(.data[["sampleID"]],
                            .data[["cellularBarcode"]],
                            sep = "_"),
            sampleID = as.factor(.data[["sampleID"]])
        ) %>%
        column_to_rownames() %>%
        .[, c("sampleID", "cellularBarcode", "nCount")]
    df
}



#' Cellular Barcode Distributions
#'
#' @author Michael Steinbaugh
#' @keywords internal
#' @noRd
#'
#' @importFrom pbapply pblapply
#'
#' @param sampleDirs Sample directories.
#'
#' @return [list].
.cellularBarcodesList <- function(sampleDirs) {
    files <- file.path(
        normalizePath(sampleDirs),
        paste(basename(sampleDirs), "barcodes.tsv", sep = "-")
    )
    names(files) <- names(sampleDirs)
    if (!all(file.exists(files))) {
        stop("Cellular barcode file missing", call. = FALSE)
    }
    message("Reading cellular barcode distributions")
    list <- pblapply(seq_along(files), function(a) {
        readFileByExtension(
            files[[a]],
            col_names = c("cellularBarcode", "nCount"),
            col_types = "ci") %>%
            mutate(cellularBarcode = gsub(
                x = .data[["cellularBarcode"]],
                pattern = "-",
                replacement = "_"))
    })
    names(list) <- names(sampleDirs)
    list
}
