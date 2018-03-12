#' Top Barcodes
#'
#' @name topBarcodes
#' @family Clustering Utilities
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#'
#' @param n Number of barcodes to return per sample.
#'
#' @return `data.frame`.
#'
#' @examples
#' load(system.file("extdata/bcb.rda", package = "bcbioSingleCell"))
#' load(system.file("extdata/seurat.rda", package = "bcbioSingleCell"))
#'
#' # bcbioSingleCell
#' topBarcodes(bcb) %>% glimpse()
#'
#' # seurat
#' topBarcodes(seurat) %>% glimpse()
NULL



# Constructors =================================================================
#' @importFrom tibble as_tibble column_to_rownames rownames_to_column
.topBarcodes <- function(object, n = 10L) {
    col <- "nUMI"
    metrics <- metrics(object)
    if (!col %in% colnames(metrics)) {
        warn(paste0("`", col, "` column missing from metrics"))
        return(NULL)
    }
    metrics %>%
        rownames_to_column() %>%
        as_tibble() %>%
        .[order(.[[col]], decreasing = TRUE), , drop = FALSE] %>%
        # Take the top rows by using slice
        dplyr::slice(1L:n) %>%
        as.data.frame() %>%
        column_to_rownames()
}



# Methods ======================================================================
#' @rdname topBarcodes
#' @export
setMethod(
    "topBarcodes",
    signature("bcbioSingleCell"),
    .topBarcodes
)




#' @rdname topBarcodes
#' @export
setMethod(
    "topBarcodes",
    signature("seurat"),
    .topBarcodes
)
