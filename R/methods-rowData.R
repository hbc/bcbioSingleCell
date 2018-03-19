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
        rowRanges <- rowRanges(x)
        assert_is_all_of(rowRanges, "GRanges")
        as(rowRanges, return)
    }
)



#' @rdname rowData
#' @export
setMethod(
    "rowData",
    signature("seurat"),
    function(x, return = c("data.frame", "DataFrame")) {
        return <- match.arg(return)
        rowRanges <- bcbio(x, "rowRanges")
        if (is(rowRanges, "GRanges")) {
            as(rowRanges, return)
        } else {
            NULL
        }
    }
)
