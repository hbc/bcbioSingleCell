#' Bind Cellular Barcodes
#'
#' @author Michael Steinbaugh
#'
#' @param list List of cellular barcodes.
#'
#' @return [tibble].
#' @noRd
.bindCellularBarcodes <- function(list) {
    lapply(seq_along(list), function(a) {
        sampleID <- names(list)[[a]]
        list[[a]] %>%
            mutate(sampleID = !!sampleID)
    }) %>%
        bind_rows() %>%
        as("tibble") %>%
        mutate(rowname = paste(.data[["sampleID"]],
                               .data[["cellularBarcode"]],
                               sep = "_"))
}
