#' Fetch Gene Expression Data
#'
#' @name fetchGeneData
#' @family Data Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#'
#' @return `matrix`.
#'
#' @examples
#' # SingleCellExperiment ====
#' genes <- rownames(bcb_small) %>% head(2L)
#' fetchGeneData(bcb_small, genes = genes) %>% glimpse()
NULL



# Methods ======================================================================
#' @rdname fetchGeneData
#' @export
setMethod(
    "fetchGeneData",
    signature("SingleCellExperiment"),
    function(object, genes) {
        assert_is_character(genes)
        counts(object) %>%
            .[genes, , drop = FALSE] %>%
            as.matrix() %>%
            t()
    }
)
