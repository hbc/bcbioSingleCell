#' Sample Metadata
#'
#' @rdname sampleMetadata
#' @name sampleMetadata
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
setMethod("sampleMetadata", "bcbioSCFiltered", .sampleMetadata)
