#' Bind Cellular Barcodes
#'
#' @rdname internal-bindCB
#' @keywords internal
#'
#' @param list List of cellular barcodes.
#'
#' @return [tibble].
.bindCB <- function(list) {
    lapply(seq_along(list), function(a) {
        sampleID <- names(list)[[a]]
        list[[a]] %>%
            mutate(sampleID = !!sampleID)
    }) %>%
        bind_rows %>%
        as("tibble") %>%
        mutate(rowname = paste(.data[["sampleID"]],
                               .data[["cellularBarcode"]],
                               sep = "_"))
}
