#' Metadata Table
#'
#' Table output for a knit report.
#'
#' @rdname metadataTable
#' @name metadataTable
#' @author Michael Steinbaugh
#'
#' @param ... Additional arguments, passed to [knitr::kable()].
#'
#' @return [kable].
NULL



# Constructors ====
.metadataTable <- function(object, ...) {
    object %>%
        sampleMetadata %>%
        mutate(fileName = NULL) %>%
        kable(caption = "Sample metadata", ...)
}



# Methods ====
#' @rdname metadataTable
#' @export
setMethod("metadataTable", "bcbioSCDataSet", .metadataTable)



#' @rdname metadataTable
#' @export
setMethod("metadataTable", "bcbioSCSubset", .metadataTable)
