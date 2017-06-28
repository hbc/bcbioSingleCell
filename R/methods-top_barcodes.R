#' Top barcodes
#'
#' @rdname top_barcodes
#'
#' @param object [bcbioSCDataSet].
#' @param n Number of barcodes to return per sample.
#'
#' @return [tibble].



#' @rdname top_barcodes
#' @usage NULL
.top_barcodes <- function(object, n = 2L) {
    metrics(object) %>%
        as_tibble %>%
        # FIXME Need to set as grouped tibble?
        top_n(n, !!sym("total_counts"))
    # FIXME Return as data.frame?
}



#' @rdname top_barcodes
#' @export
setMethod("top_barcodes", "bcbioSCDataSet", .top_barcodes)

#' @rdname top_barcodes
#' @export
setMethod("top_barcodes", "SummarizedExperiment", .top_barcodes)
