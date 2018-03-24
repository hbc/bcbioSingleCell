#' Top Markers
#'
#' @name topMarkers
#' @family Clustering Utilities
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
#' @return `tibble`.
#'
#' @examples
#' # grouped_df ====
#' topMarkers(all_markers_small, n = 1L)
NULL



# Constructors =================================================================
.topMarkers <- function(
    object,
    n = 10L,
    direction = "positive",
    coding = TRUE
) {
    .checkSanitizedMarkers(object, stop = TRUE)
    # Check to make sure `avgLogFC` column exists
    if (!"avgLogFC" %in% colnames(object)) {
        abort("`avgLogFC` column is missing")
    }
    # Check for valid direction
    directionArgs <- c("positive", "negative", "both")
    if (!direction %in% directionArgs) {
        abort(paste(
            "Valid `direction`:",
            toString(directionArgs)
        ))
    }
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
        # Arrange by P value
        .[order(.[["pvalue"]]), , drop = FALSE] %>%
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
