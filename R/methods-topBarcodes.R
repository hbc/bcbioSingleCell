#' Top Barcodes
#'
#' @rdname topBarcodes
#' @name topBarcodes
#' @family Clustering Utilities
#' @author Michael Steinbaugh
#'
#' @inheritParams AllGenerics
#' @param n Number of barcodes to return per sample.
#'
#' @return [data.frame]
#'
#' @examples
#' \dontrun{
#' data(bcb)
#' topBarcodes(bcb)
#' }
NULL



# Constructors ====
.topBarcodes <- function(object, n = 10) {
    metrics <- metrics(object) %>%
        as("tibble")
    # Check for unfiltered barcode counts in `nCount`
    if (!"nCount" %in% colnames(metrics)) {
        warning("'nCount' missing from 'metrics()'")
        return(NULL)
    }
    metrics %>%
        .[order(.[["nCount"]], decreasing = TRUE), , drop = FALSE] %>%
        # Take the top rows by using slice
        dplyr::slice(1:n) %>%
        as.data.frame() %>%
        column_to_rownames()
}



# Methods ====
#' @rdname topBarcodes
#' @export
setMethod("topBarcodes", "bcbioSingleCellANY", .topBarcodes)
