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
#' @importFrom bcbioBase checkInterestingGroups
#' @importFrom S4Vectors metadata
#' @export
setMethod(
    "interestingGroups<-",
    signature(
        object = "bcbioSingleCell",
        value = "character"),
    function(object, value) {
        sampleMetadata <- sampleMetadata(object)
        interestingGroups <- checkInterestingGroups(
            object = sampleMetadata,
            interestingGroups = value)
        metadata(object)[["interestingGroups"]] <- interestingGroups
        validObject(object)
        object
    })



#' @rdname interestingGroups
#' @importFrom bcbioBase checkInterestingGroups
#' @importFrom S4Vectors metadata
#' @export
setMethod(
    "interestingGroups<-",
    signature(
        object = "seurat",
        value = "character"),
    function(object, value) {
        sampleMetadata <- sampleMetadata(object)
        interestingGroups <- checkInterestingGroups(
            object = sampleMetadata,
            interestingGroups = value)
        bcbio(object, "interestingGroups") <- interestingGroups
        validObject(object)
        object
    })
