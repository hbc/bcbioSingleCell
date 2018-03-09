# TODO Need to add `rowRanges()` support for seurat



#' Row Data
#'
#' @name rowData
#' @author Michael Steinbaugh
#'
#' @importFrom SummarizedExperiment rowData
#'
#' @inheritParams general
#'
#' @return Return data as `data.frame`, or `DataFrame`.
#'
#' @examples
#' load(system.file("extdata/bcb.rda", package = "bcbioSingleCell"))
#' load(system.file("extdata/seurat.rda", package = "bcbioSingleCell"))
#'
#' # bcbioSingleCell ====
#' rowData(bcb) %>% glimpse()
#'
#' # seurat ====
#' rowData(pbmc_small) %>% glimpse()
#' rowData(seurat) %>% glimpse()
NULL



# Methods ======================================================================
#' @rdname rowData
#' @export
setMethod(
    "rowData",
    signature("bcbioSingleCell"),
    function(x, return = c("data.frame", "DataFrame")) {
        return <- match.arg(return)
        data <- mcols(rowRanges(x))
        rownames(data) <- names(rowRanges(x))
        as(data, return)
    }
)



#' @rdname rowData
#' @export
setMethod(
    "rowData",
    signature("seurat"),
    function(x, return = c("data.frame", "DataFrame")) {
        return <- match.arg(return)
        # Catch `rowData` or `annotable`
        # TODO Error on stashed `annotable` in future update
        match <- match(
            x = c("rowData", "annotable"),
            table = names(bcbio(x))
        ) %>%
            na.omit()
        if (!length(match)) {
            return(NULL)
        } else {
            match <- match[[1L]]
        }
        data <- bcbio(x)[[match]]
        as(data, return)
    }
)
