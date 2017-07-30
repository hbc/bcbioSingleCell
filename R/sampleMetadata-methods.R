#' Sample Metadata
#'
#' @rdname sampleMetadata
#' @name sampleMetadata
#' @author Michael Steinbaugh
#'
#' @return [tibble].
NULL



# Constructors ====
.sampleMetadata <- function(object) {
    metadata(object)[["sampleMetadata"]]
}



# Methods ====
#' @rdname sampleMetadata
#' @export
setMethod("sampleMetadata", "bcbioSCDataSet", .sampleMetadata)



#' @rdname sampleMetadata
#' @export
setMethod("sampleMetadata", "bcbioSCSubset", .sampleMetadata)
