#' Sample Metadata
#'
#' @rdname sampleMetadata
#' @name sampleMetadata
#' @author Michael Steinbaugh
#'
#' @return [data.frame].
NULL



# Constructors ====
.sampleMetadata <- function(object) {
    metadata(object)[["sampleMetadata"]] %>%
        as.data.frame
}



# Methods ====
#' @rdname sampleMetadata
#' @export
setMethod("sampleMetadata", "bcbioSCDataSet", .sampleMetadata)



#' @rdname sampleMetadata
#' @export
setMethod("sampleMetadata", "bcbioSCSubset", .sampleMetadata)
