#' Sample Metadata
#'
#' @rdname sample_metadata
#' @author Michael Steinbaugh
#'
#' @return [tibble].



#' @rdname sample_metadata
#' @usage NULL
.sample_metadata <- function(object) {
    metadata(object)[["sample_metadata"]]
}



#' @rdname sample_metadata
#' @export
setMethod("sample_metadata", "bcbioSCDataSet", .sample_metadata)



#' @rdname sample_metadata
#' @export
setMethod("sample_metadata", "bcbioSCSubset", .sample_metadata)
