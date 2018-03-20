#' Fetch Gene Expression Data
#'
#' @name fetchGeneData
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#' @param genes Gene identifiers (matrix rownames).
#'
#' @return `matrix`.
#'
#' @examples
#' # seurat ====
#' genes <- counts(seurat_small) %>% rownames() %>% head()
#' fetchGeneData(seurat_small, genes = genes) %>% glimpse()
NULL



# Methods ======================================================================
#' @rdname fetchGeneData
#' @importFrom Seurat FetchData
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
