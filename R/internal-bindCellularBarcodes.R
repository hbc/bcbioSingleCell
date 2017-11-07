#' Bind Cellular Barcodes
#'
#' @author Michael Steinbaugh
#' @keywords internal
#' @noRd
#'
#' @importFrom dplyr bind_rows mutate mutate_if
#' @importFrom parallel mclapply
#' @importFrom tibble as_tibble
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
        as_tibble() %>%
        mutate(
            cellID = paste(.data[["sampleID"]],
                           .data[["cellularBarcode"]],
                           sep = "_")
        ) %>%
        mutate_if(is.character, as.factor)
}
