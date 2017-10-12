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
#' @export
setMethod(
    "annotable",
    signature("bcbioSingleCellANY"),
    function(object) {
        metadata(object)[["annotable"]]
    })
