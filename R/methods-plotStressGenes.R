#' Plot Stress Genes
#'
#' @rdname plotStressGenes
#' @name plotStressGenes
#' @family Clustering Utilities
#' @author Michael Steinbaugh
#' @keywords internal
#'
#' @inheritParams AllGenerics
#'
#' @return [ggplot].
NULL



# Constructors ====
# Immediate early genes
.ieg <- function() {
}



# Integrated stress response
.isr <- function() {
}



# Methods ====
#' @rdname plotStressGenes
#' @export
setMethod(
    "plotStressGenes",
    signature("ANY"),
    function(object) {
    stop("Draft function")
})
