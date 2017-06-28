#' Counts accessor
#'
#' @rdname counts
#'
#' @author Michael Steinbaugh
#'
#' @param object [bcbioSCDataSet].
#' @param format Matrix format.
#'
#' @return Counts matrix.



#' @rdname counts
#' @usage NULL
.counts <- function(object, format = "dgCMatrix") {
    supported_formats <- c("dgCMatrix", "dgTMatrix", "matrix")
    if (!format %in% supported_formats) {
        stop(paste("Supported formats", toString(supported_formats)))
    }
    message(format)
    assay(object) %>% as(format)
}



#' @rdname counts
#' @export
setMethod("counts", "bcbioSCDataSet", .counts)

#' @rdname counts
#' @export
setMethod("counts", "SummarizedExperiment", .counts)
