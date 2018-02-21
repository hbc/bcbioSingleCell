#' Detect Organism
#'
#' @rdname detectOrganism
#' @name detectOrganism
#' @author Michael Steinbaugh
#'
#' @importFrom basejump detectOrganism
#'
#' @inheritParams general
#'
#' @examples
#' load(system.file("extdata/bcb.rda", package = "bcbioSingleCell"))
#' load(system.file("extdata/seurat.rda", package = "bcbioSingleCell"))
#'
#' # bcbioSingleCell
#' detectOrganism(bcb)
#'
#' # seurat
#' detectOrganism(seurat)
NULL



# Methods ======================================================================
#' @rdname detectOrganism
#' @export
setMethod(
    "detectOrganism",
    signature("bcbioSingleCell"),
    function(object) {
        organism <- metadata(object)[["organism"]]
        if (is.character(organism)) {
            detectOrganism(organism)
        } else {
            warn(paste(
                "`organism` is not defined in metadata() slot.",
                "Attempting to match by the first gene identifier row",
                "in the counts matrix instead."
            ))
            assay(object) %>%
                rownames() %>%
                .[[1L]] %>%
                detectOrganism()
        }
    }
)


#' @rdname detectOrganism
#' @export
setMethod(
    "detectOrganism",
    signature("seurat"),
    function(object) {
        organism <- slot(object, "misc") %>%
            .[["bcbio"]] %>%
            .[["organism"]]
        if (is.character(organism)) {
            detectOrganism(organism)
        } else {
            warn(paste(
                "`organism` is not defined in `bcbio() slot`.",
                "Attempting to match by the first gene identifier row",
                "in the counts matrix instead."
            ))
            slot(object, "data") %>%
                rownames() %>%
                .[[1L]] %>%
                detectOrganism()
        }
    }
)
