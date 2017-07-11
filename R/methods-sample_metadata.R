#' Sample metadata
#'
#' @rdname sample_metadata
#'
#' @author Michael Steinbaugh
#'
#' @param object Primary object.
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
setMethod("sample_metadata", "SCSubset", .sample_metadata)
