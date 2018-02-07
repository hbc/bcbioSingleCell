#' Update Object
#'
#' @note SingleCellExperiment migration started in v0.0.31.
#'
#' @rdname updateObject
#' @name updateObject
#'
#' @author Michael Steinbaugh
#'
#' @inheritParams AllGenerics
#'
#' @examples
#' loadRemoteData("http://bcbiosinglecell.seq.cloud/v0.0.26/bcb.rda")
#'
#' # bcbioSingleCell prior to v0.0.31 extended from SummarizedExperiment instead
#' of SingleCellExperiment
#' metadata(bcb)[["version"]]
#' error <- tryCatch(
#'     validObject(bcb),
#'     error = function(e) e
#' )
#' print(error)
#' bcbSCE <- updateObject(bcb)
#' validObject(bcbSCE)
#' is(bcbSCE, "SingleCellExperiment")
#' bcbSCE
NULL



# Methods ======================================================================
#' @rdname updateObject
#' @export
setMethod(
    "updateObject",
    signature("bcbioSingleCell"),
    function(object) {
        version <- metadata(object)[["version"]]
        if (version < "0.0.31") {
            se <- as(object, "SummarizedExperiment")
            sce <- as(se, "SingleCellExperiment")
            bcbio <- bcbio(object)
            new("bcbioSingleCell", sce, bcbio = bcbio)
        } else {
            object
        }
    })
