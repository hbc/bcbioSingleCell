#' Top barcodes
#'
#' @param run [bcbioSCDataSet].
#' @param n Number of barcodes to return per sample.
#'
#' @export
top_barcodes <- function(run, n = 2) {
    run$metrics %>%
        top_n(n, !!sym("total_counts"))
}
