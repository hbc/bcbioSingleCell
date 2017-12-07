# FIXME Add Seurat support



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
#' # bcbioSingleCell
#' bcb <- examples[["bcb"]]
#' topBarcodes(bcb) %>% glimpse()
#'
#' # seurat
#' seurat <- examples[["seurat"]]
NULL



# Constructors ====
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



# Methods ====
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
