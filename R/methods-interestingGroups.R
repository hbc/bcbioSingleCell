#' Interesting Groups
#'
#' @name interestingGroups
#' @family Data Functions
#' @author Michael Steinbaugh
#'
#' @importFrom bcbioBase interestingGroups interestingGroups<-
#'
#' @inheritParams general
#'
#' @return Character vector.
#'
#' @examples
#' load(system.file("extdata/bcb_small.rda", package = "bcbioSingleCell"))
#' load(system.file("extdata/seurat_small.rda", package = "bcbioSingleCell"))
#'
#' # bcbioSingleCell ====
#' interestingGroups(bcb_small)
#'
#' # Assignment method support
#' interestingGroups(bcb_small) <- "sampleID"
#' interestingGroups(bcb_small)
#'
#' # seurat ====
#' interestingGroups(seurat_small)
#'
#' # Assignment method support
#' interestingGroups(seurat_small) <- "sampleID"
#' interestingGroups(seurat_small)
NULL



# Constructors =================================================================
.interestingGroups <- function(object) {
    validObject(object)
    x <- metadata(object)[["interestingGroups"]]
    if (is.null(x)) {
        x <- "sampleName"
    }
    x
}



`.interestingGroups<-` <- function(object, value) {
    assertFormalInterestingGroups(
        x = sampleData(object),
        interestingGroups = value
    )
    metadata(object)[["interestingGroups"]] <- value
    validObject(object)
    object
}



# Methods ======================================================================
#' @rdname interestingGroups
#' @export
setMethod(
    "interestingGroups",
    signature("bcbioSingleCell"),
    .interestingGroups
)



#' @rdname interestingGroups
#' @export
setMethod(
    "interestingGroups",
    signature("seurat"),
    .interestingGroups
)



# Assignment methods ===========================================================
#' @rdname interestingGroups
#' @export
setMethod(
    "interestingGroups<-",
    signature(
        object = "bcbioSingleCell",
        value = "character"
    ),
    `.interestingGroups<-`
)



#' @rdname interestingGroups
#' @export
setMethod(
    "interestingGroups<-",
    signature(
        object = "seurat",
        value = "character"
    ),
    `.interestingGroups<-`
)
