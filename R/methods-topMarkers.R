#' Top Markers
#'
#' @rdname topMarkers
#' @name topMarkers
#' @family Clustering Utilities
#' @author Michael Steinbaugh
#'
#' @inheritParams AllGenerics
#'
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
#' @importFrom dplyr arrange filter slice
#' @importFrom rlang !! sym
.topMarkers <- function(
    object,
    n = 10,
    coding = TRUE) {
    .validMarkers(object)
    if (isTRUE(coding)) {
        object <- object %>%
            filter(.data[["biotype"]] == "protein_coding") %>%
            # Remove additional genes annotated as "predicted" in description
            filter(!grepl(x = .data[["description"]],
                          pattern = "^predicted\\s"))
    }
    object %>%
        # Use only the positive markers
        filter(.data[["avgDiff"]] > 0) %>%
        # Arrange by P value
        arrange(!!sym("pvalue")) %>%
        # Take the top rows by using slice
        slice(1:n)
}



# Methods ====
#' @rdname topMarkers
#' @export
setMethod(
    "topMarkers",
    signature("grouped_df"),
    .topMarkers)
