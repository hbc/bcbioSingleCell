#' Row Data
#'
#' Metadata describing the genes. These annotations are obtained from
#' AnnotationHub and ensembldb using [basejump::ensembl()].
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
#' # bcbioSingleCell ====
#' rowData(bcb_small) %>% glimpse()
#'
#' # seurat ====
#' rowData(pbmc_small) %>% glimpse()
#' rowData(seurat_small) %>% glimpse()
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
