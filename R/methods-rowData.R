#' Row Data
#'
#' @name rowData
#' @author Michael Steinbaugh
#'
#' @importFrom SummarizedExperiment rowData
#'
#' @inheritParams general
#'
#' @return `DataFrame`, `data.frame`, or `GRanges`.
#'
#' @examples
#' load(system.file("extdata/bcb.rda", package = "bcbioSingleCell"))
#' load(system.file("extdata/seurat.rda", package = "bcbioSingleCell"))
#'
#' # bcbioSingleCell
#' rowData(bcb) %>% glimpse()
#'
#' # seurat
#' rowData(seurat) %>% glimpse()
NULL



# Constructors =================================================================
.rowData <- function(x, return = c("DataFrame", "data.frame", "AsIs")) {
    return <- match.arg(return)
    data <- slot(x, "elementMetadata")
    if (return != "AsIs") {
        data <- as(data, return)
    }
    names <- slot(x, "NAMES")
    if (has_dims(data)) {
        rownames(data) <- names
    } else if (has_names(data)) {
        names(data) <- names
    }
    data
}



# Methods ======================================================================
#' @rdname rowData
#' @export
setMethod(
    "rowData",
    signature("bcbioSingleCell"),
    .rowData
)



#' @rdname rowData
#' @export
setMethod(
    "rowData",
    signature("seurat"),
    function(x, return = c("DataFrame", "data.frame", "AsIs")) {
        return <- match.arg(return)
        data <- bcbio(x, "rowData")
        if (return != "AsIs") {
            data <- as(data, return)
        }
        data
    }
)
