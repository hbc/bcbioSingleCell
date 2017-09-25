#' Read Cellular Barcode File
#'
#' @author Michael Steinbaugh
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
        mutate(cellularBarcode = str_replace_all(
            .data[["cellularBarcode"]], "-", "_"))
}
