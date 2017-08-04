#' Interesting Groups
#'
#' @rdname interestingGroups
#' @name interestingGroups
#'
#' @return Character vector.
NULL



# Constructors ====
.interestingGroups <- function(object) {
    metadata(object)[["interestingGroups"]]
}



# Methods ====
#' @rdname interestingGroups
#' @export
setMethod("interestingGroups", "bcbioSCDataSet", .interestingGroups)



#' @rdname interestingGroups
#' @export
setMethod("interestingGroups", "bcbioSCSubset", .interestingGroups)
