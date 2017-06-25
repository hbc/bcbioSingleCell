#' Sample barcode metrics
#'
#' @rdname metrics
#'
#' @param object [bcbioSCDataSet].
#'
#' @return [data.frame].
#' @export
setMethod("metrics", "bcbioSCDataSet", function(object) {
    colData(object) %>% as.data.frame
})
