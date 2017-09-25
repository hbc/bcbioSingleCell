#' Known Markers Detected
#'
#' @rdname knownMarkersDetected
#' @name knownMarkersDetected
#' @family Clustering Utilities
#' @author Michael Steinbaugh
#'
#' @param object [data.frame]: [Seurat::FindAllMarkers()] return.
#' @param knownMarkers Known markers [data.frame] imported by [readMarkers()] or
#'   pulled from internal [cellCycleData] object.
#' @param show Show [kable].
#'
#' @return [tibble].
NULL



# Constructors ====
.knownMarkersDetected <- function(object, knownMarkers, show = FALSE) {
    markers <- .groupMarkers(object) %>%
        left_join(knownMarkers[, c("cell", "ensgene")],
                  by = "ensgene") %>%
        dplyr::select(c("cell", "ensgene", "symbol", "cluster"),
                    everything()) %>%
        dplyr::filter(!is.na(.data[["cell"]])) %>%
        .[order(.[["cell"]],
                .[["symbol"]],
                .[["pvalue"]],
                .[["avgDiff"]],
                decreasing = c(FALSE, FALSE, FALSE, TRUE)), ] %>%
        group_by(!!sym("cell"))
    if (isTRUE(show)) {
        .markersKable(markers, caption = "Known markers detected")
    }
    markers
}



# Methods ====
#' @rdname knownMarkersDetected
#' @export
setMethod(
    "knownMarkersDetected",
    signature(object = "data.frame", knownMarkers = "tbl_df"),
    .knownMarkersDetected)
