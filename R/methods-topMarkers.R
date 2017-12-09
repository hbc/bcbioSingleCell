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
#' @param direction Whether to include only `positive`, `negative`, or `both`
#'   directions of association per cluster. Defaults to `positive`.
#' @param coding Only include protein coding genes.
#'
#' @seealso
#' - [dplyr::slice()]
#' - [dplyr::top_n()].
#'
#' @return [tibble].
#' @export
NULL



# Constructors =================================================================
#' @importFrom dplyr arrange filter slice
#' @importFrom rlang !! sym
.topMarkers <- function(
    object,
    n = 10,
    direction = "positive",
    coding = TRUE) {
    .checkSanitizedMarkers(object, stop = TRUE)
    # Check to make sure `avgLogFC` column exists
    if (!"avgLogFC" %in% colnames(object)) {
        stop("'avgLogFC' column is required")
    }
    # Check for valid direction
    directionArgs <- c("positive", "negative", "both")
    if (!direction %in% directionArgs) {
        stop(paste(
            "Valid 'direction':",
            toString(directionArgs)
        ), call. = FALSE)
    }
    if (isTRUE(coding)) {
        object <- object %>%
            filter(.data[["biotype"]] == "protein_coding") %>%
            # Remove additional genes annotated as "predicted" in description
            filter(!grepl(x = .data[["description"]],
                          pattern = "^predicted\\s"))
    }
    # Subset to positive or negative correlation, if desired ("direction")
    # Note that `avgDiff` has been renamed to `avgLogFC` in Seurat v2.1
    if (direction == "positive") {
        object <- filter(object, .data[["avgLogFC"]] > 0)
    } else if (direction == "negative") {
        object <- filter(object, .data[["avgLogFC"]] < 0)
    }
    object %>%
        # Arrange by P value
        arrange(!!sym("pvalue")) %>%
        # Take the top rows by using slice
        slice(1:n)
}



# Methods ======================================================================
#' @rdname topMarkers
#' @export
setMethod(
    "topMarkers",
    signature("grouped_df"),
    .topMarkers)
