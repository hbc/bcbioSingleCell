#' Top Barcodes
#'
#' @rdname topBarcodes
#' @name topBarcodes
#' @family Clustering Utilities
#' @author Michael Steinbaugh
#'
#' @inheritParams AllGenerics
#'
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
#' @importFrom dplyr slice
#' @importFrom tibble as_tibble column_to_rownames rownames_to_column
.topBarcodes <- function(object, n = 10) {
    metrics <- metrics(object) %>%
        rownames_to_column() %>%
        as_tibble()
    # Check for unfiltered barcode counts in `nCount`
    if (!"nCount" %in% colnames(metrics)) {
        warning("'nCount' missing from 'metrics()'", call. = FALSE)
        return(NULL)
    }
    metrics %>%
        .[order(.[["nCount"]], decreasing = TRUE), , drop = FALSE] %>%
        # Take the top rows by using slice
        slice(1:n) %>%
        as.data.frame() %>%
        column_to_rownames()
}



# Methods ====
#' @rdname topBarcodes
#' @export
setMethod(
    "topBarcodes",
    signature("bcbioSingleCell"),
    .topBarcodes)
