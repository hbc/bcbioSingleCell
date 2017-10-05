#' Bind Cellular Barcodes
#'
#' @author Michael Steinbaugh
#' @keywords internal
#' @noRd
#'
#' @param list List of cellular barcodes.
#'
#' @return [tibble].
.bindCellularBarcodes <- function(list) {
    mclapply(seq_along(list), function(a) {
        sampleID <- names(list)[[a]]
        list[[a]] %>%
            mutate(sampleID = !!sampleID)
    }) %>%
        bind_rows() %>%
        as("tibble") %>%
        mutate(cellID = paste(.data[["sampleID"]],
                              .data[["cellularBarcode"]],
                              sep = "_"))
}
