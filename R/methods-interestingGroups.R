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
#' @examples
#' # bcbioSingleCell ====
#' interestingGroups(bcb_small)
#' interestingGroups(bcb_small) <- "sampleName"
#' interestingGroups(bcb_small)
NULL



# Methods ======================================================================
#' @rdname interestingGroups
#' @export
setMethod(
    "interestingGroups<-",
    signature(
        object = "SingleCellExperiment",
        value = "character"
    ),
    getMethod(
        "interestingGroups<-",
        signature(
            object = "SummarizedExperiment",
            value = "character"
        )
    )
)
