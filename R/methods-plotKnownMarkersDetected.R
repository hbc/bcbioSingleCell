#' Plot Known Markers Detected
#'
#' @rdname plotKnownMarkersDetected
#' @name plotKnownMarkersDetected
#' @family Clustering Utilities
#' @author Michael Steinbaugh
#'
#' @inherit plotMarkers
#'
#' @param knownMarkersDetected [knownMarkersDetected()] return `grouped_df`,
#'   grouped by `cluster` ident.
#'
#' @examples
#' # seurat ====
#' plotKnownMarkersDetected(
#'     object = seurat_small,
#'     knownMarkersDetected = known_markers_small
#' )
NULL



# Constructors =================================================================
#' @importFrom basejump markdownHeader
.plotKnownMarkersDetected <- function(
    object,
    knownMarkersDetected,
    tsneColor = scale_color_viridis(),
    violinFill = scale_fill_viridis(discrete = TRUE),
    dotColor = scale_color_gradient(
        low = "lightgray",
        high = "purple"
    ),
    dark = TRUE,
    pointsAsNumbers = FALSE,
    headerLevel = 2L
) {
    # Passthrough: tsneColor, violinFill, dotColor, dark, pointsAsNumbers
    assertIsAHeaderLevel(headerLevel)
    assert_has_rows(knownMarkersDetected)
    assert_is_subset("cellType", colnames(knownMarkersDetected))

    cellTypes <- knownMarkersDetected %>%
        pull("cellType") %>%
        na.omit() %>%
        unique()
    assert_is_non_empty(cellTypes)

    list <- pblapply(cellTypes, function(cellType) {
        genes <- knownMarkersDetected %>%
            .[.[["cellType"]] == cellType, , drop = FALSE] %>%
            pull("geneName") %>%
            na.omit() %>%
            unique()
        assert_is_non_empty(genes)

        markdownHeader(
            cellType,
            level = headerLevel,
            tabset = TRUE,
            asis = TRUE
        )
        subheaderLevel <- headerLevel + 1L

        plotMarkers(
            object,
            genes = genes,
            tsneColor = tsneColor,
            violinFill = violinFill,
            dotColor = dotColor,
            dark = dark,
            pointsAsNumbers = pointsAsNumbers,
            headerLevel = subheaderLevel,
            title = cellType
        )
    })

    invisible(list)
}



# Methods ======================================================================
#' @rdname plotKnownMarkersDetected
#' @export
setMethod(
    "plotKnownMarkersDetected",
    signature(
        object = "seurat",
        knownMarkersDetected = "grouped_df"
    ),
    .plotKnownMarkersDetected
)
