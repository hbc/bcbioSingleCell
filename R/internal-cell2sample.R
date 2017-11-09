#' Define Cell to Sample Mappings
#'
#' This function extracts `sampleID` from the `cellID` column using grep
#' matching.
#'
#' @importFrom dplyr bind_rows filter mutate
#' @importFrom magrittr set_colnames
#' @importFrom parallel mclapply
#' @importFrom stringr str_match
#'
#' @author Michael Steinbaugh
#' @keywords internal
#' @noRd
.cell2sample <- function(cells, samples) {
    list <- mclapply(seq_along(samples), function(a) {
        pattern <- paste0("^(", samples[[a]], barcodePattern)
        match <- str_match(cells, pattern = pattern) %>%
            as.data.frame() %>%
            set_colnames(
                c("cellID",
                  "sampleID",
                  "cellularBarcode",
                  # Trailing number is for matching cellranger
                  "trailingNumber")
            ) %>%
            filter(!is.na(.data[["sampleID"]])) %>%
            mutate(
                cellularBarcode = paste0(
                    .data[["cellularBarcode"]],
                    na.omit(.data[["trailingNumber"]])
                ),
                trailingNumber = NULL
            )
        match
    })
    cell2sample <- bind_rows(list)
    if (!identical(length(cells), nrow(cell2sample))) {
        stop("Failed to correctly match sample IDs to cellular barcodes",
             call. = FALSE)
    }
    cell2sample
}
