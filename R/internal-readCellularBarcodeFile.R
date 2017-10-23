#' Read Cellular Barcode File
#'
#' @author Michael Steinbaugh
#'
#' @importFrom basejump readFileByExtension
#' @importFrom dplyr mutate
#'
#' @param file Cellular barcode TSV file.
#'
#' @return [tibble].
#' @noRd
.readCellularBarcodeFile <- function(file) {
    readFileByExtension(
        file,
        col_names = c("cellularBarcode", "nCount"),
        col_types = "ci") %>%
        mutate(cellularBarcode = gsub(
            x = .data[["cellularBarcode"]],
            pattern = "-",
            replacement = "_"))
}
