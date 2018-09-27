#' Counts per Cellular Barcode
#'
#' @author Michael Steinbaugh
#' @keywords internal
#' @noRd
#'
#' @param list `list`. Cellular barcodes per sample.
#'
#' @return `integer`.
.nCount <- function(list) {
    assert_is_list(list)
    assert_has_names(list)
    assert_is_integer(list[[1L]])
    assert_has_names(list[[1L]])
    if (length(list) > 1L) {
        # Prefix with `sampleID` when handling multiple samples.
        mapply(
            sampleID = names(list),
            x = list,
            FUN = function(x, sampleID) {
                # Prefix the barcode with `sampleID`.
                names(x) <- paste(sampleID, names(x), sep = "_")
                x
            },
            USE.NAMES = TRUE,
            SIMPLIFY = TRUE
        )
    } else {
        list[[1L]]
    }
}
