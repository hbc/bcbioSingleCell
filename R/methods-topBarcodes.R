#' Top Barcodes
#'
#' @rdname topBarcodes
#' @name topBarcodes
#' @family Clustering Utilities
#' @author Michael Steinbaugh
#'
#' @param n Number of barcodes to return per sample.
#'
#' @return [tibble].
NULL



# Constructors ====
.topBarcodes <- function(object, n = 2L) {
    metrics(object) %>%
        as("tibble") %>%
        # Use slice here instead?
        top_n(n, !!sym("total_counts"))
}



# Methods ====
#' @rdname topBarcodes
#' @export
setMethod("topBarcodes", "bcbioSCDataSet", .topBarcodes)



#' @rdname topBarcodes
#' @export
setMethod("topBarcodes", "bcbioSCFiltered", .topBarcodes)
