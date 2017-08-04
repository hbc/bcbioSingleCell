#' Top Barcodes
#'
#' @rdname topBarcodes
#' @name topBarcodes
#'
#' @param n Number of barcodes to return per sample.
#'
#' @return [tibble].
NULL



# Constructors ====
.topBarcodes <- function(object, n = 2L) {
    metrics(object) %>%
        as("tibble") %>%
        top_n(n, !!sym("total_counts"))
}



# Methods ====
#' @rdname topBarcodes
#' @export
setMethod("topBarcodes", "bcbioSCDataSet", .topBarcodes)



#' @rdname topBarcodes
#' @export
setMethod("topBarcodes", "bcbioSCSubset", .topBarcodes)
