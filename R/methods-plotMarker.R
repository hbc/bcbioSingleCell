#' Plot Cell-Type Gene Marker(s)
#'
#' @description
#' Plot gene expression per cell in multiple formats:
#'
#' 1. [plotMarkerTSNE()]: t-SNE gene expression plot.
#' 2. [plotDot()]: Dot plot.
#' 3. [plotViolin()]: Violin plot.
#'
#' @name plotMarker
#' @family Clustering Functions
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @inheritParams general
#'
#' @return Show graphical output. Invisibly return `ggplot` `list`.
NULL



# Constructors =================================================================
# Strip everything except the x-axis text labels
.minimalAxis <- function() {
    theme(
        axis.line = element_blank(),
        # axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        panel.grid = element_blank(),
        title = element_blank()
    )
}



# Methods ======================================================================
#' @rdname plotMarker
#' @inheritParams fetchTSNEExpressionData
#' @inheritParams plotTSNE
#' @param expression Calculation to apply on the aggregate marker expression.
#' @export
#' @examples
#' plotMarkerTSNE(seurat_small, genes = "COL1A1")
#' plotMarkerTSNE(seurat_small, genes = "COL1A1", dark = FALSE, grid = FALSE)
#'
#' # Mitochondrial genes
#' mito <- grep("^MT\\.", rownames(counts(seurat_small)), value = TRUE)
#' print(sort(mito))
#' plotMarkerTSNE(seurat_small, genes = mito, title = "mitochondrial")
setMethod(
    "plotMarkerTSNE",
    signature("seurat"),
    function(
        object,
        genes,
        expression = c("mean", "median", "sum"),
        color = "auto",
        pointsAsNumbers = FALSE,
        pointSize = 0.5,
        pointAlpha = 0.8,
        label = TRUE,
        labelSize = 6L,
        dark = TRUE,
        grid = TRUE,
        legend = FALSE,
        aspectRatio = 1L,
        title = TRUE
    ) {
        assert_is_character(genes)
        assert_is_subset(genes, rownames(object))
        expression <- match.arg(expression)
        assert_is_a_bool(pointsAsNumbers)
        assert_is_a_number(pointSize)
        assert_is_a_number(pointAlpha)
        assert_is_a_bool(label)
        assert_is_a_number(labelSize)
        assert_is_a_bool(dark)
        assert_is_a_bool(legend)
        assert_is_any_of(title, c("character", "logical", "NULL"))
        if (is.character(title)) {
            assert_is_a_string(title)
        }

        data <- fetchTSNEExpressionData(object, genes = genes)
        requiredCols <- c(
            "centerX",
            "centerY",
            "mean",
            "median",
            "ident",
            "sum",
            "tSNE1",
            "tSNE2"
        )
        assert_is_subset(requiredCols, colnames(data))

        p <- ggplot(
            data = data,
            mapping = aes_string(
                x = "tSNE1",
                y = "tSNE2",
                color = expression
            )
        )

        # Titles
        subtitle <- NULL
        if (isTRUE(title)) {
            if (is_a_string(genes)) {
                title <- genes
            } else {
                title <- NULL
                subtitle <- genes
                # Limit to the first 5 markers
                if (length(subtitle) > 5L) {
                    subtitle <- c(subtitle[1L:5L], "...")
                }
                subtitle <- toString(subtitle)
            }
        } else if (identical(title, FALSE)) {
            title <- NULL
        }
        p <- p + labs(title = title, subtitle = subtitle)

        # Customize legend
        if (isTRUE(legend)) {
            # Make the guide longer than normal, to improve appearance of values
            # containing a decimal point
            p <- p +
                guides(
                    color = guide_colorbar(
                        barwidth = 20L,
                        barheight = 1L,
                        direction = "horizontal"
                    )
                )
        } else {
            p <- p + guides(color = "none")
        }

        if (isTRUE(pointsAsNumbers)) {
            p <- p +
                geom_text(
                    mapping = aes_string(
                        x = "tSNE1",
                        y = "tSNE2",
                        label = "ident",
                        color = expression
                    ),
                    alpha = pointAlpha,
                    size = pointSize
                )
        } else {
            p <- p +
                geom_point(
                    alpha = pointAlpha,
                    size = pointSize
                )
        }

        if (isTRUE(label)) {
            if (isTRUE(dark)) {
                labelColor <- "white"
            } else {
                labelColor <- "black"
            }
            p <- p +
                geom_text(
                    mapping = aes_string(
                        x = "centerX",
                        y = "centerY",
                        label = "ident"
                    ),
                    color = labelColor,
                    size = labelSize,
                    fontface = "bold"
                )
        }

        # Color palette
        if (isTRUE(dark)) {
            p <- p +
                theme_midnight(
                    aspect_ratio = aspectRatio,
                    grid = grid
                )
            if (color == "auto") {
                color <- scale_color_viridis(
                    option = "plasma",
                    discrete = FALSE
                )
            }
        } else {
            p <- p +
                theme_paperwhite(
                    aspect_ratio = aspectRatio,
                    grid = grid
                )
            if (color == "auto") {
                color <- scale_color_gradient(
                    low = "gray90",
                    high = "black"
                )
            }
        }

        if (is(color, "ScaleContinuous")) {
            p <- p + color
        }

        p
    }
)



#' @rdname plotMarker
#' @param gene Gene identifier. Must intersect with [rownames()].
#' @export
#' @examples
#' # Individual gene
#' plotMarker(seurat_small, gene = "COL1A2", dark = TRUE)
#' plotMarker(seurat_small, gene = "COL1A2", dark = FALSE)
setMethod(
    "plotMarker",
    signature("seurat"),
    function(
        object,
        gene,
        dark = TRUE,
        return = c("grid", "list"),
        ...
    ) {
        assert_is_a_string(gene)
        gene <- make.names(gene)
        assert_is_a_bool(dark)
        return <- match.arg(return)

        # Plots ================================================================
        if (isTRUE(dark)) {
            violinFill <- "white"
        } else {
            violinFill <- "black"
        }

        tsne <- plotMarkerTSNE(
            object,
            genes = gene,
            expression = "sum",
            dark = dark,
            ...
        )

        dot <- plotDot(
            object,
            genes = gene,
            dark = dark
        )

        violin <- plotViolin(
            object,
            genes = gene,
            scale = "width",
            fill = violinFill,
            dark = dark,
            return = "list"
        )
        # Get the ggplot object from the list return
        violin <- violin[[1L]]

        # Return ===============================================================
        if (return == "grid") {
            violin <- violin +
                .minimalAxis()
            dot <- dot +
                coord_flip() +
                .minimalAxis()
            p <- plot_grid(
                tsne,
                dot,
                violin,
                labels = NULL,
                ncol = 1L,
                nrow = 3L,
                rel_heights = c(1L, 0.1, 0.15)
            )
            if (isTRUE(dark)) {
                p <- p + theme(
                    plot.background = element_rect(color = NA, fill = "black")
                )
            } else {
                p <- p + theme(
                    plot.background = element_rect(color = NA, fill = "white")
                )
            }
            p
        } else if (return == "list") {
            list(
                "tsne" = tsne,
                "dot" = dot,
                "violin" = violin
            )
        }
    }
)



#' @rdname plotMarker
#' @export
#' @examples
#' # Multiple genes
#' top <- topMarkers(all_markers_small, n = 1L)
#' genes <- pull(top, "rowname")
#' plotMarkers(seurat_small, genes = genes)
setMethod(
    "plotMarkers",
    signature("seurat"),
    function(
        object,
        genes,
        headerLevel = 2L,
        ...
    ) {
        assert_is_subset(genes, rownames(object))
        assertIsAHeaderLevel(headerLevel)

        list <- lapply(genes, function(gene) {
            markdownHeader(gene, level = headerLevel, asis = TRUE)
            p <- plotMarker(
                object = object,
                gene = gene,
                ...
            )
            show(p)
            invisible(p)
        })

        invisible(list)
    }
)



#' @rdname plotMarker
#' @note The number of markers to plot is determined by the output of the
#' [topMarkers()] function. If you want to reduce the number of genes to plot,
#' simply reassign first using that function. If necessary, we can add support
#' for the number of genes to plot here in a future update.
#' @param topMarkers `grouped_df` grouped by "`cluster`", returned by
#'   [topMarkers()].
#' @export
#' @examples
#' plotTopMarkers(
#'     object = seurat_small,
#'     topMarkers = topMarkers(all_markers_small, n = 1L)
#' )
setMethod(
    "plotTopMarkers",
    signature("seurat"),
    function(
        object,
        topMarkers,
        headerLevel = 2L,
        ...
    ) {
        validObject(object)
        stopifnot(is(topMarkers, "grouped_df"))
        stopifnot(.isSanitizedMarkers(topMarkers))
        assertIsAHeaderLevel(headerLevel)

        clusters <- levels(topMarkers[["cluster"]])
        list <- pblapply(clusters, function(cluster) {
            genes <- topMarkers %>%
                .[.[["cluster"]] == cluster, , drop = FALSE] %>%
                pull("rowname")
            if (!length(genes)) {
                return(invisible())
            }
            if (length(genes) > 10L) {
                warning("Maximum of 10 genes per cluster is recommended")
            }
            markdownHeader(
                paste("Cluster", cluster),
                level = headerLevel,
                tabset = TRUE,
                asis = TRUE
            )
            subheaderLevel <- headerLevel + 1L
            plotMarkers(
                object = object,
                genes = genes,
                headerLevel = subheaderLevel,
                ...
            )
        })

        invisible(list)
    }
)



#' @rdname plotMarker
#' @param markers `grouped_df` of known marker genes.
#' @export
#' @examples
#' # Let's plot the first known marker, as a quick example
#' plotKnownMarkersDetected(
#'     object = seurat_small,
#'     markers = known_markers_small[1L, , drop = FALSE]
#' )
setMethod(
    "plotKnownMarkersDetected",
    signature("seurat"),
    function(
        object,
        markers,
        headerLevel = 2L,
        ...
    ) {
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
                object = cellType,
                level = headerLevel,
                tabset = TRUE,
                asis = TRUE
            )
            subheaderLevel <- headerLevel + 1L

            plotMarkers(
                object = object,
                genes = genes,
                headerLevel = subheaderLevel,
                ...
            )
        })

        invisible(list)
    }
)
