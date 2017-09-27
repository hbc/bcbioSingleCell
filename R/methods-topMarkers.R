#' Top Markers
#'
#' @rdname topMarkers
#' @name topMarkers
#' @family Clustering Utilities
#' @author Michael Steinbaugh
#'
#' @inheritParams AllGenerics
#' @param n Number of genes per cluster.
#' @param coding Only include protein coding genes.
#'
#' @seealso
#' - [dplyr::slice()]
#' - [dplyr::top_n()].
#'
#' @return [tibble].
#' @export
NULL



# Constructors ====
.topMarkers <- function(
    object,
    n = 4L,
    coding = FALSE,
    show = FALSE) {
    .validMarkers(object)
    if (isTRUE(coding)) {
        object <- dplyr::filter(
            object, .data[["biotype"]] == "protein_coding")
    }
    object %>%
        # Use only the positive markers
        dplyr::filter(.data[["avgDiff"]] > 0L) %>%
        # Arrange by P value
        arrange(!!sym("pvalue")) %>%
        # Take the top rows by using slice
        dplyr::slice(1:n)
}



# Methods ====
#' @rdname topMarkers
#' @export
setMethod("topMarkers", "grouped_df", .topMarkers)
