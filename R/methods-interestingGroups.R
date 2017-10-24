#' Interesting Groups
#'
#' @rdname interestingGroups
#' @name interestingGroups
#' @author Michael Steinbaugh
#'
#' @inheritParams AllGenerics
#'
#' @return Character vector.
NULL



# Methods ====
#' @rdname interestingGroups
#' @export
setMethod(
    "interestingGroups",
    signature("bcbioSingleCell"),
    function(object) {
        metadata(object)[["interestingGroups"]]
    })
