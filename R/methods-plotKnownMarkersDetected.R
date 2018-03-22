#' Plot Known Markers Detected
#'
#' @rdname plotKnownMarkersDetected
#' @name plotKnownMarkersDetected
#' @family Clustering Utilities
#' @author Michael Steinbaugh
#'
#' @inherit plotMarkers
#'
#' @param knownMarkersDetected `grouped_df`.
#'
#' @examples
#' load(system.file(
#'     "extdata/knownMarkersDetected.rda",
#'     package = "bcbioSingleCell"
#' ))
#' load(system.file("extdata/seurat.rda", package = "bcbioSingleCell"))
#'
#' # seurat ====
#' # Let's plot the first 2 markers, as a quick example
#' plotKnownMarkersDetected(
#'     seurat,
#'     knownMarkersDetected = knownMarkersDetected[seq_len(2L), ]
#' )
NULL



# Constructors =================================================================
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
    if (!nrow(knownMarkersDetected)) {
        return(NULL)
    }
    cellTypes <- knownMarkersDetected[, "cellType", drop = TRUE] %>%
        na.omit()
        unique()
    if (is.null(cellTypes)) return(NULL)
    list <- pblapply(seq_along(cellTypes), function(a) {
        cellType <- cellTypes[[a]]
        genes <- knownMarkersDetected %>%
            .[.[["cellType"]] == cellType, "geneName", drop = TRUE] %>%
            na.omit() %>%
            unique()
        if (is.null(genes)) {
            return(NULL)
        }
        if (!is.null(headerLevel)) {
            markdownHeader(
                cellType,
                level = headerLevel,
                tabset = TRUE,
                asis = TRUE
            )
            subheaderLevel <- headerLevel + 1L
        } else {
            subheaderLevel <- NULL
        }
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
