#' Top Markers
#'
#' @name topMarkers
#' @family Clustering Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#' @param n Number of genes per cluster.
#' @param direction Whether to include "`positive`", "`negative`", or "`both`"
#'   directions of association per cluster.
#' @param coding Only include protein coding genes.
#'
#' @seealso
#' - [dplyr::slice()].
#' - [dplyr::top_n()].
#'
#' @return `grouped_df`.
#'
#' @examples
#' # grouped_df ====
#' topMarkers(
#'     object = all_markers_small,
#'     n = 2L,
#'     direction = "positive"
#' )
NULL



# Constructors =================================================================
.topMarkers <- function(
    object,
    n = 10L,
    direction = c("positive", "negative", "both"),
    coding = TRUE
) {
    stopifnot(.isSanitizedMarkers(object))
    assert_is_subset(c("avgLogFC", "padj"), colnames(object))
    assertIsAnImplicitInteger(n)
    direction <- match.arg(direction)
    assert_is_a_bool(coding)

    if (isTRUE(coding)) {
        object <- object %>%
            .[.[["geneBiotype"]] == "protein_coding", , drop = FALSE] %>%
            # Remove additional genes annotated as "predicted" in description
            .[!grepl("^predicted\\s", .[["description"]]), , drop = FALSE]
    }

    # Subset to positive or negative correlation, if desired ("direction")
    # Note that `avgDiff` has been renamed to `avgLogFC` in Seurat v2.1
    if (direction == "positive") {
        object <- object[object[["avgLogFC"]] > 0L, , drop = FALSE]
    } else if (direction == "negative") {
        object <- object[object[["avgLogFC"]] < 0L, , drop = FALSE]
    }

    object %>%
        # Arrange by adjusted P value
        arrange(!!sym("padj")) %>%
        # Take the top rows by using slice
        dplyr::slice(1L:n)
}



# Methods ======================================================================
#' @rdname topMarkers
#' @export
setMethod(
    "topMarkers",
    signature("grouped_df"),
    .topMarkers
)
