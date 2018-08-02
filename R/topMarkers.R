#' Top Markers
#'
#' @family Clustering Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#' @param data `grouped_df`. `Seurat::FindAllMarkers()` return, sanitized with
#'   our `sanitizeSeuratMarkers()` function.
#' @param n `scalar integer`. Number of genes per cluster.
#' @param direction `string`. Whether to include "`positive`", "`negative`", or
#'   "`both`" directions of association per cluster.
#' @param coding `boolean`. Include only protein coding genes.
#'
#' @seealso
#' - [dplyr::slice()].
#' - [dplyr::top_n()].
#'
#' @return `grouped_df`.
#' @export
#'
#' @examples
#' # grouped_df ====
#' x <- topMarkers(
#'     data = all_markers_small,
#'     n = 2L,
#'     direction = "positive"
#' )
#' head(x)
topMarkers <- function(
    data,
    n = 10L,
    direction = c("positive", "negative", "both"),
    coding = FALSE
) {
    stopifnot(.isSanitizedMarkers(data))
    assert_is_subset(c("avgLogFC", "padj"), colnames(data))
    assertIsAnImplicitInteger(n)
    direction <- match.arg(direction)
    assert_is_a_bool(coding)

    if (isTRUE(coding)) {
        assert_is_subset("geneBiotype", colnames(data))
        data <- data %>%
            .[.[["geneBiotype"]] == "protein_coding", , drop = FALSE] %>%
            # Remove genes annotated as "predicted"
            .[!grepl("^predicted\\s", .[["description"]]), , drop = FALSE]
    }

    # Subset to positive or negative correlation, if desired ("direction")
    # Note that `avgDiff` has been renamed to `avgLogFC` in Seurat v2.1
    if (direction == "positive") {
        data <- data[data[["avgLogFC"]] > 0L, , drop = FALSE]
    } else if (direction == "negative") {
        data <- data[data[["avgLogFC"]] < 0L, , drop = FALSE]
    }

    data %>%
        # Arrange by adjusted P value
        arrange(!!sym("padj")) %>%
        # Take the top rows by using slice
        dplyr::slice(1L:n)
}
