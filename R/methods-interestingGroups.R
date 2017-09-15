#' Interesting Groups
#'
#' @rdname interestingGroups
#' @name interestingGroups
#' @author Michael Steinbaugh
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
setMethod("interestingGroups", "bcbioSCFiltered", .interestingGroups)
