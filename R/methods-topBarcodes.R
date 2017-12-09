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
#' load(system.file(
#'     file.path("extdata", "bcb.rda"),
#'     package = "bcbioSingleCell"))
#' load(system.file(
#'     file.path("extdata", "seurat.rda"),
#'     package = "bcbioSingleCell"))
#'
#' # bcbioSingleCell
#' topBarcodes(bcb) %>% glimpse()
#'
#' # seurat
#' topBarcodes(seurat) %>% glimpse()
NULL



# Constructors =================================================================
#' @importFrom dplyr slice
#' @importFrom tibble as_tibble column_to_rownames rownames_to_column
.topBarcodes <- function(object, n = 10) {
    col <- "nUMI"
    metrics <- metrics(object)
    if (!col %in% colnames(metrics)) {
        warning(paste0(
            "'", col, "' missing from 'metrics()'"
        ), call. = FALSE)
        return(NULL)
    }
    metrics %>%
        rownames_to_column() %>%
        as_tibble() %>%
        .[order(.[[col]], decreasing = TRUE), , drop = FALSE] %>%
        # Take the top rows by using slice
        slice(1:n) %>%
        as.data.frame() %>%
        column_to_rownames()
}



# Methods ======================================================================
#' @rdname topBarcodes
#' @export
setMethod(
    "topBarcodes",
    signature("bcbioSingleCell"),
    .topBarcodes)




#' @rdname topBarcodes
#' @export
setMethod(
    "topBarcodes",
    signature("seurat"),
    .topBarcodes)
