#' Ensembl Annotations
#'
#' @rdname annotable
#' @name annotable
#' @author Michael Steinbaugh
#'
#' @inheritParams AllGenerics
#'
#' @return [data.frame]
NULL



# Methods ====
#' @rdname annotable
#' @importFrom S4Vectors metadata
#' @export
setMethod(
    "annotable",
    signature("bcbioSingleCell"),
    function(object) {
        metadata(object)[["annotable"]]
    })
