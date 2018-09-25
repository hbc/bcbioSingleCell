#' Bind Cellular Barcodes
#'
#' @author Michael Steinbaugh
#' @keywords internal
#' @noRd
#'
#' @param list `list`. Cellular barcodes per sample.
#'
#' @return `data.frame`.
.bindCellularBarcodes <- function(list) {
    assert_is_list(list)
    assert_has_names(list)
    mapply(
        sampleID = names(list),
        x = list,
        FUN = function(x, sampleID) {
            assert_are_identical(colnames(x), c("barcode", "n"))
            # Add the sampleID as a column
            x[["sampleID"]] <- as.factor(sampleID)
            x
        },
        USE.NAMES = TRUE,
        SIMPLIFY = FALSE
    ) %>%
        bind_rows() %>%
        mutate(
            rowname = paste(!!sym("sampleID"), !!sym("barcode"), sep = "_")
        ) %>%
        as("DataFrame")
}
