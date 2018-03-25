#' Plot Known Markers Detected
#'
#' @name plotKnownMarkersDetected
#' @family Clustering Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams plotMarkers
#' @inheritParams general
#' @param knownMarkersDetected `grouped_df`.
#'
#' @return Show graphical output. Invisibly return `ggplot` plotlist.
#'
#' @examples
#' load(system.file(
#'     "extdata/known_markers_small.rda",
#'     package = "bcbioSingleCell"
#' ))
#' load(system.file("extdata/seurat_small.rda", package = "bcbioSingleCell"))
#'
#' # seurat ====
#' # Let's plot the first 2 markers, as a quick example
#' plotKnownMarkersDetected(
#'     object = seurat_small,
#'     knownMarkersDetected = known_markers_small[seq_len(2L), ]
#' )
NULL



# Constructors =================================================================
.plotKnownMarkersDetected <- function(
    object,
    knownMarkersDetected,
    tsneColor = scale_color_viridis(discrete = FALSE),
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
    assert_has_rows(knownMarkersDetected)
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
            object = object,
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
