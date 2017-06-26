#' Top barcodes
#'
#' @rdname top_barcodes
#'
#' @param object [bcbioSCDataSet].
#' @param n Number of barcodes to return per sample.
#' @param ... Additional parameters.
#'
#' @return [data.frame].
#' @export
setMethod("top_barcodes", "bcbioSCDataSet", function(object, n = 2L) {
    metrics(object) %>%
        as_tibble %>%
        # FIXME Need to set as grouped tibble?
        top_n(n, !!sym("total_counts"))
    # FIXME Return as data.frame?
})
