#' Top Barcodes
#'
#' @rdname top_barcodes
#' @author Michael Steinbaugh
#'
#' @param n Number of barcodes to return per sample.
#'
#' @return [tibble].



#' @rdname top_barcodes
#' @usage NULL
.top_barcodes <- function(object, n = 2L) {
    metrics(object) %>%
        as("tibble") %>%
        top_n(n, !!sym("total_counts"))
}



#' @rdname top_barcodes
#' @export
setMethod("top_barcodes", "bcbioSCDataSet", .top_barcodes)

#' @rdname top_barcodes
#' @export
setMethod("top_barcodes", "bcbioSCSubset", .top_barcodes)
