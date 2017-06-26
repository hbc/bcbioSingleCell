#' Counts accessor
#'
#' @rdname counts
#'
#' @author Michael Steinbaugh
#'
#' @param object Object.
#' @param ... Additional parameters.
#'
#' @return Sparse counts (`dgCMatrix`).
#' @export
setMethod("counts", "bcbioSCDataSet", function(object) {
    assay(object)
})
