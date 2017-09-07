#' Read Cellular Barcode File
#'
#' @author Michael Steinbaugh
#' @keywords internal
#'
#' @param file Cellular barcode TSV file.
#'
#' @return [tibble].
.readCBFile <- function(file) {
    readFileByExtension(
        file,
        col_names = c("cellularBarcode", "nCount"),
        col_types = "ci") %>%
        mutate(cellularBarcode = str_replace_all(
            .data[["cellularBarcode"]], "-", "_"))
}
