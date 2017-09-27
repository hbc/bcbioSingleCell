#' Metadata
#'
#' @rdname metadata
#' @name metadata
#' @author Michael Steinbaugh
#' @keywords internal
NULL



# Methods ====
#' @rdname metadata
#' @export
setMethod("metadata", "bcbioSingleCell", function(x) {
    x@metadata
})
