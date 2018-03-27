#' Gene to Symbol Mappings
#'
#' @name gene2symbol
#' @family Data Functions
#' @author Michael Steinbaugh
#'
#' @importFrom basejump gene2symbol
#'
#' @inheritParams general
#'
#' @examples
#' # bcbioRNASeq ====
#' gene2symbol(bcb_small) %>% glimpse()
#'
#' # seurat ====
#' gene2symbol(seurat_small) %>% glimpse()
NULL



# Methods ======================================================================
#' @rdname gene2symbol
#' @export
setMethod(
    "gene2symbol",
    signature("seurat"),
    function(object) {
        rowData <- rowData(object)
        assert_is_non_empty(rowData)
        cols <- c("geneID", "geneName")
        assert_is_subset(cols, colnames(rowData))
        rowData[, cols]
    }
)
