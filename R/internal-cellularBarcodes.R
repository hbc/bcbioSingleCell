#' Bind Cellular Barcodes
#'
#' @author Michael Steinbaugh
#' @keywords internal
#' @noRd
#'
#' @importFrom dplyr bind_rows group_by mutate mutate_if
#' @importFrom parallel mclapply
#' @importFrom rlang !!
#' @importFrom tibble as_tibble
#'
#' @param list List of cellular barcodes.
#'
#' @return [tibble] grouped by `sampleID`.
.bindCellularBarcodes <- function(list) {
    list <- mclapply(seq_along(list), function(a) {
        # Add the sampleID as a column
        sampleID <- names(list)[[a]]
        list[[a]] %>%
            mutate(sampleID = !!sampleID)
    })
    tbl <- list %>%
        bind_rows() %>%
        as_tibble() %>%
        mutate(
            cellID = paste(.data[["sampleID"]],
                           .data[["cellularBarcode"]],
                           sep = "_")
        ) %>%
        .[, c("sampleID", "cellularBarcode", "cellID", "nCount")] %>%
        mutate_if(is.character, as.factor) %>%
        group_by(.data[["sampleID"]])
    tbl
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
    files <- file.path(paste(
        basename(sampleDirs),
        "barcodes.tsv",
        sep = "-"
    ))
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
