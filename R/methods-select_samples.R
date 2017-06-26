# TODO Switch to AnnotationDbi::select method?

#' Select samples
#'
#' @rdname select_samples
#'
#' @param object Object.
#' @param samples Character vector of desired samples.
#' @param ... Additional parameters.
#'
#' @return [bcbioSCDataSet].
#' @export
setMethod("select_samples", "bcbioSCDataSet", function(object, samples) {
    samples <- c("6189_sample", "6191_sample")
    pattern <- paste0(
        "^(",
        paste(samples, collapse = "|"),
        ")")
    assay(object) %>%
        .[, str_detect(colnames(sparse_counts), pattern)]
})
