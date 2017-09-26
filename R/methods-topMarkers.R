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
#' @param show Show [kable].
#'
#' @seealso
#' - [dplyr::slice()]
#' - [dplyr::top_n()].
#'
#' @return [tibble].
#' @export
NULL



# Constructors ====
# Also called by `knownMarkersDetected()`
.markersKable <- function(object, caption = NULL) {
    object %>%
        .[, c("cluster",
              "symbol",
              "ensgene",
              "pvalue",
              "description",
              "biotype")] %>%
        # Format the P values into consistent scientific notation
        mutate(pvalue = format(.data[["pvalue"]],
                               digits = 3L,
                               scientific = TRUE)) %>%
        # Remove the source information from description
        mutate(description = str_replace(
            .data[["description"]], " \\[.+$", "")) %>%
        # Truncate the description to 50 characters
        mutate(description = str_trunc(.data[["description"]], 50L)) %>%
        kable(caption = caption) %>%
        show()
}



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
    object <- object %>%
        # Use only the positive markers
        dplyr::filter(.data[["avgDiff"]] > 0L) %>%
        # Arrange by P value
        arrange(!!sym("pvalue")) %>%
        # Take the top rows by using slice
        dplyr::slice(1:n)
    if (isTRUE(show)) {
        .markersKable(object, caption = paste("Top", n, "markers per cluster"))
    }
    object
}



# Methods ====
#' @rdname topMarkers
#' @export
setMethod("topMarkers", "grouped_df", .topMarkers)
