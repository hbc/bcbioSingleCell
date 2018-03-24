#' Fetch Gene Expression Data
#'
#' @name fetchGeneData
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#'
#' @return `matrix`.
#'
#' @examples
#' # seurat ====
#' genes <- head(rownames(pbmc_small))
#' fetchGeneData(pbmc_small, genes = genes) %>% glimpse()
NULL



# Methods ======================================================================
#' @rdname fetchGeneData
#' @export
setMethod(
    "fetchGeneData",
    signature("seurat"),
    function(object, genes) {
        assert_is_character(genes)
        data <- Seurat::FetchData(object, vars.all = genes)
        assert_are_identical(
            x = as.character(genes),
            y = as.character(colnames(data))
        )
        data
    }
)
