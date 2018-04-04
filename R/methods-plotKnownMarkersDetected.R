#' Plot Known Markers Detected
#'
#' @name plotKnownMarkersDetected
#' @family Clustering Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams plotMarkers
#' @inheritParams general
#' @param markers `grouped_df` of known marker genes.
#'
#' @return Show graphical output. Invisibly return `ggplot` `list`.
#'
#' @examples
#' # seurat ====
#' # Let's plot the first 2 markers, as a quick example
#' plotKnownMarkersDetected(
#'     object = seurat_small,
#'     markers = known_markers_small[seq_len(2L), ]
#' )
NULL



# Methods ======================================================================
#' @rdname plotKnownMarkersDetected
#' @export
setMethod(
    "plotKnownMarkersDetected",
    signature("seurat"),
    function(
        object,
        markers,
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
        assert_has_rows(markers)
        assertIsAHeaderLevel(headerLevel)
        stopifnot(is(markers, "grouped_df"))
        assert_has_rows(markers)
        assert_is_subset("cellType", colnames(markers))

        cellTypes <- markers %>%
            pull("cellType") %>%
            na.omit() %>%
            unique()
        assert_is_non_empty(cellTypes)

        list <- pblapply(cellTypes, function(cellType) {
            genes <- markers %>%
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
)
