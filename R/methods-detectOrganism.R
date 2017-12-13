#' Detect Organism
#'
#' @rdname detectOrganism
#' @name detectOrganism
#' @author Michael Steinbaugh
#'
#' @importFrom basejump detectOrganism
#' @inherit basejump::detectOrganism
#'
#' @examples
#' load(system.file(
#'     file.path("extdata", "bcb.rda"),
#'     package = "bcbioSingleCell"))
#' load(system.file(
#'     file.path("extdata", "seurat.rda"),
#'     package = "bcbioSingleCell"))
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
            warning(paste(
                "'organism' is not slotted in metadata(object)",
                "Attempting to match by the first gene identifier row",
                "in the counts matrix instead."
            ), call. = FALSE)
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
            warning(paste(
                "'organism' is not slotted in 'bcbio(object)'.",
                "Attempting to match by the first gene identifier row",
                "in the counts matrix instead."
                ), call. = FALSE)
            slot(object, "data") %>%
                rownames() %>%
                .[[1]] %>%
                detectOrganism()
        }
    }
)
