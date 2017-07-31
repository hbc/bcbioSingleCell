#' Bind Cellular Barcodes
#'
#' @rdname bindCB-internal
#' @family Cellular Barcode Utilities
#' @author Michael Steinbaugh
#' @keywords internal
#'
#' @param list List of cellular barcodes.
#'
#' @return [tibble].
.bindCB <- function(list) {
    lapply(seq_along(list), function(a) {
        sampleID <- names(list)[[a]] %>% camel
        list[[a]] %>%
            mutate(sampleID = !!sampleID)
    }) %>%
        as("tibble") %>%
        bind_rows %>%
        mutate(rowname = paste(.data[["sampleID"]],
                               .data[["cellularBarcode"]],
                               sep = "_"))
}
