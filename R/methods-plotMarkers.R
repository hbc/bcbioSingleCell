# TODO Add `plotMarkerUMAP()`, dimensional reduction option to `plotMarkers()`



#' Plot Cell-Type Gene Marker(s)
#'
#' @description
#' Plot gene expression per cell in multiple formats:
#'
#' 1. [plotMarkerTSNE()]: t-SNE gene expression plot.
#' 2. [plotDot()]: Dot plot.
#' 3. [plotViolin()]: Violin plot.
#'
#' @section plotTopMarkers:
#' The number of markers to plot is determined by the output of the
#' [topMarkers()] function. If you want to reduce the number of genes to plot,
#' simply reassign first using that function. If necessary, we can add support
#' for the number of genes to plot here in a future update.
#'
#' @name plotMarkers
#' @family Clustering Functions
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @inheritParams fetchTSNEExpressionData
#' @inheritParams plotTSNE
#' @inheritParams general
#' @param expression Calculation to apply on the aggregate marker expression.
#' @param markers `grouped_df` of marker genes.
#'   - [plotTopMarkers()]: must be grouped by "`cluster`".
#'   - [plotKnownMarkersDetected()]: must be grouped by "`cellType`".
#'
#' @return Show graphical output. Invisibly return `ggplot` `list`.
#'
#' @examples
#' object <- seurat_small
#' genes <- "COL1A1"
#'
#' # plotMarkerTSNE ====
#' plotMarkerTSNE(object, genes = genes)
#' plotMarkerTSNE(object, genes = genes, dark = FALSE, grid = FALSE)
#'
#' # Mitochondrial genes
#' mito <- grep("^MT-", rownames(object), value = TRUE)
#' print(mito)
#' plotMarkerTSNE(object, genes = mito, title = "mito")
#'
#' # plotMarkers ====
#' plotMarkers(object, genes = genes)
#'
#' # plotTopMarkers ====
#' markers <- topMarkers(all_markers_small, n = 1)
#' markers
#' markers <- head(markers, n = 1)
#' plotTopMarkers(object, markers = markers)
#'
#' # plotKnownMarkersDetected ====
#' markers <- head(known_markers_small, n = 1)
#' markers
#' plotKnownMarkersDetected(object, markers = markers)
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



.plotMarker <- function(
    object,
    dimRed = c("tsne", "umap"),
    gene,
    dark = TRUE,
    ...
) {
    dimRed <- match.arg(dimRed)

    if (isTRUE(dark)) {
        violinFill <- "white"
    } else {
        violinFill <- "black"
    }

    # Dimensional reduction plot
    if (dimRed == "tsne") {
        dr <- plotMarkerTSNE(
            object,
            genes = gene,
            expression = "sum",
            dark = dark,
            ...
        )
    } else if (dimRed == "umap") {
        # FIXME Add UMAP support
        stop("UMAP not supported yet")
    }

    # Dot plot
    dot <- plotDot(
        object,
        genes = gene,
        dark = dark
    ) +
        coord_flip() +
        .minimalAxis()

    # Violin plot
    violin <- plotViolin(
        object,
        genes = gene,
        scale = "width",
        fill = violinFill,
        dark = dark,
        return = "list"
    ) %>%
        # Get the ggplot object from the list return
        .[[1L]] +
        .minimalAxis()

    plotlist <- list(
        dr = dr,
        dot = dot,
        violin = violin
    )

    # Return ===============================================================
    p <- plot_grid(
        plotlist = plotlist,
        labels = NULL,
        ncol = 1L,
        nrow = 3L,
        rel_heights = c(1L, 0.1, 0.15)
    )
    if (isTRUE(dark)) {
        p <- p + theme(
            plot.background = element_rect(
                color = NA,
                fill = "black"
            )
        )
    } else {
        p <- p + theme(
            plot.background = element_rect(
                color = NA,
                fill = "white"
            )
        )
    }
    p
}



# Methods ======================================================================
#' @rdname plotMarkers
#' @export
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



#' @rdname plotMarkers
#' @export
setMethod(
    "plotMarkers",
    signature("seurat"),
    function(
        object,
        genes,
        dark = TRUE,
        headerLevel = 2L,
        ...
    ) {
        assert_is_subset(genes, rownames(object))
        assert_is_a_bool(dark)
        assertIsAHeaderLevel(headerLevel)
        list <- lapply(genes, function(gene) {
            p <- .plotMarker(object = object, gene = gene, ...)
            markdownHeader(gene, level = headerLevel, asis = TRUE)
            show(p)
            invisible(p)
        })
        names(list) <- genes
        invisible(list)
    }
)



#' @rdname plotMarkers
#' @export
setMethod(
    "plotTopMarkers",
    signature("seurat"),
    function(
        object,
        markers,
        headerLevel = 2L,
        ...
    ) {
        validObject(object)
        stopifnot(is(markers, "grouped_df"))
        stopifnot(.isSanitizedMarkers(markers))
        assertIsAHeaderLevel(headerLevel)

        clusters <- levels(markers[["cluster"]])
        list <- pblapply(clusters, function(cluster) {
            genes <- markers %>%
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



#' @rdname plotMarkers
#' @export
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
