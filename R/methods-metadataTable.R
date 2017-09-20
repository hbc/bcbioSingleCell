#' Metadata Table
#'
#' Table output for a knit report.
#'
#' @rdname metadataTable
#' @name metadataTable
#' @author Michael Steinbaugh
#'
#' @inheritParams AllGenerics
#' @param ... Additional arguments, passed to [knitr::kable()].
#'
#' @return [kable].
NULL



# Constructors ====
.metadataTable <- function(object, ...) {
    sampleMetadata(object) %>%
        mutate(fileName = NULL) %>%
        # Put sample name first and sort
        dplyr::select("sampleName", everything()) %>%
        arrange(!!sym("sampleName")) %>%
        kable(caption = "Sample metadata", ...)
}



# Methods ====
#' @rdname metadataTable
#' @export
setMethod("metadataTable", "bcbioSCDataSet", .metadataTable)



#' @rdname metadataTable
#' @export
setMethod("metadataTable", "bcbioSCFiltered", .metadataTable)
