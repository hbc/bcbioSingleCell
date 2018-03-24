#' Counts Accessor
#'
#' @name counts
#' @author Michael Steinbaugh
#'
#' @importFrom BiocGenerics counts
#'
#' @inheritParams general
#' @param normalized Normalized (`TRUE`) or raw (`FALSE`) counts.
#'
#' @return Counts matrix (typically a `dgCMatrix`).
#'
#' @examples
#' load(system.file("extdata/bcb_small.rda", package = "bcbioSingleCell"))
#'
#' # bcbioSingleCell ====
#' counts(bcb_small) %>% glimpse()
#'
#' # seurat ====
#' counts(pbmc_small, normalized = FALSE) %>% glimpse()
#' counts(pbmc_small, normalized = TRUE) %>% glimpse()
NULL



# Methods ======================================================================
#' @rdname counts
#' @importFrom basejump makeNames
#' @export
setMethod(
    "counts",
    signature("bcbioSingleCell"),
    function(object) {
        assay(object)
    }
)



#' @rdname counts
#' @export
setMethod(
    "counts",
    signature("seurat"),
    function(object, normalized = FALSE) {
        assert_is_a_bool(normalized)
        # seurat also stashes scaled counts in `scale.data`
        if (identical(normalized, FALSE)) {
            slot(object, "raw.data")
        } else if (identical(normalized, TRUE)) {
            slot(object, "data")
        }
    }
)
