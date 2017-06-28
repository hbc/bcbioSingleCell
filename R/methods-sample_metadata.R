#' Sample metadata
#'
#' @rdname sample_metadata
#'
#' @author Michael Steinbaugh
#'
#' @param object [bcbioSCDataSet].
#' @param unique_names Unique sample names.
#'
#' @return [data.frame].



#' @rdname sample_metadata
#' @usage NULL
.sample_metadata <- function(object, unique_names = FALSE) {
    meta <- metadata(object)[["sample_metadata"]] %>% as.data.frame
    if (isTRUE(unique_names)) {
        # Ensure unique sample names
        if (any(duplicated(meta[["sample_name"]]))) {
            meta[["sample_name"]] <-
                paste0(meta[["sample_name"]],
                       " (",
                       meta[["file_name"]],
                       ")")
        }
        # Double check
        if (any(duplicated(meta[["sample_name"]]))) {
            stop("Automatic generation of unique sample names failed")
        }
    }
    meta
}



#' @rdname sample_metadata
#' @export
setMethod("sample_metadata", "bcbioSCDataSet", .sample_metadata)



#' @rdname sample_metadata
#' @export
setMethod("sample_metadata", "SummarizedExperiment", .sample_metadata)
