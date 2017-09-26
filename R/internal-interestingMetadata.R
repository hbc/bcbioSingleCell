#' Interesting Metadata
#'
#' @author Michael Steinbaugh
#'
#' @inheritParams AllGenerics
#'
#' @return [data.frame].
#' @noRd
.interestingMetadata <- function(object) {
    sampleMetadata(object) %>%
        .[, unique(c(metaPriorityCols,
                     interestingGroups(object))),
          drop = FALSE]
}
