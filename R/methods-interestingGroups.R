#' Interesting Groups
#'
#' @rdname interestingGroups
#' @name interestingGroups
#' @author Michael Steinbaugh
#'
#' @importFrom bcbioBase interestingGroups interestingGroups<-
#'
#' @inheritParams general
#'
#' @return Character vector.
#'
#' @examples
#' load(system.file("extdata/bcb.rda", package = "bcbioSingleCell"))
#' load(system.file("extdata/seurat.rda", package = "bcbioSingleCell"))
#'
#' # bcbioSingleCell
#' interestingGroups(bcb)
#'
#' # Assignment method support
#' interestingGroups(bcb) <- "sampleID"
#' interestingGroups(bcb)
#'
#' # seurat
#' interestingGroups(seurat)
#'
#' # Assignment method support
#' interestingGroups(seurat) <- "sampleID"
#' interestingGroups(seurat)
NULL



# Methods ======================================================================
#' @rdname interestingGroups
#' @export
setMethod(
    "interestingGroups",
    signature("bcbioSingleCell"),
    function(object) {
        metadata(object)[["interestingGroups"]]
    })



#' @rdname interestingGroups
#' @export
setMethod(
    "interestingGroups",
    signature("seurat"),
    function(object) {
        interestingGroups <- bcbio(object)[["interestingGroups"]]
        if (is.null(interestingGroups)) {
            interestingGroups <- "sampleName"
        }
        interestingGroups
    })



# Assignment methods ===========================================================
#' @rdname interestingGroups
#' @export
setMethod(
    "interestingGroups<-",
    signature(
        object = "bcbioSingleCell",
        value = "character"),
    function(object, value) {
        assertFormalInterestingGroups(
            x = sampleData(object),
            interestingGroups = value)
        metadata(object)[["interestingGroups"]] <- value
        validObject(object)
        object
    })



#' @rdname interestingGroups
#' @export
setMethod(
    "interestingGroups<-",
    signature(
        object = "seurat",
        value = "character"),
    function(object, value) {
        assertFormalInterestingGroups(
            x = sampleData(object),
            interestingGroups = value)
        bcbio(object, "interestingGroups") <- value
        validObject(object)
        object
    })
