#' Plot Cell-Type-Specific Gene Markers
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
#' @name plotMarker
#' @family Clustering Functions
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @inheritParams general
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
#' # Markdown
#' plotMarker(object, genes = genes)
#'
#' # t-SNE
#' plotMarkerTSNE(
#'     object = object,
#'     genes = genes,
#'     dark = TRUE,
#'     grid = TRUE
#' )
#' plotMarkerTSNE(
#'     object = object,
#'     genes = genes,
#'     dark = FALSE,
#'     grid = FALSE
#' )
#'
#' # Aggregate marker expression (e.g. mitochondrial genes)
#' mito_genes <- grep("^MT-", rownames(object), value = TRUE)
#' print(mito_genes)
#' plotMarkerTSNE(
#'     object = object,
#'     genes = mito_genes,
#'     expression = "mean",
#'     title = "mito genes"
#' )
#'
#' # UMAP
#' plotMarkerUMAP(
#'     object = object,
#'     genes = genes,
#'     dark = TRUE,
#'     grid = TRUE
#' )
#' plotMarkerUMAP(
#'     object = object,
#'     genes = genes,
#'     dark = FALSE,
#'     grid = FALSE
#' )
#'
#' # Top markers
#' markers <- topMarkers(all_markers_small, n = 1)
#' markers
#' markers <- head(markers, n = 1)
#' plotTopMarkers(object, markers = markers)
#'
#' # Known markers detected
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
    gene,
    dimRed = c("tsne", "umap"),
    dark = TRUE,
    grid = TRUE
) {
    dimRed <- match.arg(dimRed)

    if (isTRUE(dark)) {
        fill <- "black"
        violinFill <- "white"
    } else {
        fill <- "white"
        violinFill <- "black"
    }

    # Dimensional reduction plot
    if (dimRed == "tsne") {
        plotMarkerDimRed <- plotMarkerTSNE
    } else if (dimRed == "umap") {
        plotMarkerDimRed <- plotMarkerUMAP
    }
    dr <- plotMarkerDimRed(
        object,
        genes = gene,
        dark = dark,
        grid = grid,
        legend = FALSE
    )

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
    plot_grid(
        plotlist = plotlist,
        labels = NULL,
        ncol = 1L,
        nrow = 3L,
        rel_heights = c(1L, 0.1, 0.15)
    ) +
        theme(
            plot.background = element_rect(
                color = NA,
                fill = fill
            )
        )
}



# Dimensionality reduction (t-SNE, UMAP) plot constructor
.plotMarkerDimRed <- function(
    object,
    genes,
    dimRed = c("tsne", "umap"),
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
    dimRed <- match.arg(dimRed)
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

    # Fetch dimensional reduction coordinates
    if (dimRed == "tsne") {
        fun <- fetchTSNEExpressionData
        dimCols <- c("tSNE1", "tSNE2")
    } else if (dimRed == "umap") {
        fun <- fetchUMAPExpressionData
        dimCols <- c("umap1", "umap2")
    }
    data <- fun(object, genes = genes)

    requiredCols <- c(
        "centerX",
        "centerY",
        "mean",
        "median",
        "ident",
        "sum",
        dimCols
    )
    assert_is_subset(requiredCols, colnames(data))

    p <- ggplot(
        data = data,
        mapping = aes_string(
            x = dimCols[[1L]],
            y = dimCols[[2L]],
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
                    x = dimCols[[1L]],
                    y = dimCols[[2L]],
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
            color <- scale_colour_viridis(option = "plasma")
        }
    } else {
        p <- p +
            theme_paperwhite(
                aspect_ratio = aspectRatio,
                grid = grid
            )
        if (color == "auto") {
            color <- scale_colour_viridis(
                option = "plasma", begin = 1L, end = 0L
            )
        }
    }

    if (is(color, "ScaleContinuous")) {
        p <- p + color
    }

    p
}



# Methods ======================================================================
#' @rdname plotMarker
#' @export
setMethod(
    "plotMarker",
    signature("seurat"),
    function(
        object,
        genes,
        dimRed = c("tsne", "umap"),
        dark = TRUE,
        grid = TRUE,
        headerLevel = 2L
    ) {
        assert_is_subset(genes, rownames(object))
        dimRed <- match.arg(dimRed)
        assert_is_a_bool(dark)
        assert_is_a_bool(grid)
        assertIsAHeaderLevel(headerLevel)
        list <- lapply(genes, function(gene) {
            p <- .plotMarker(
                object = object,
                gene = gene,
                dimRed = dimRed,
                dark = dark,
                grid = grid
            )
            markdownHeader(gene, level = headerLevel, asis = TRUE)
            show(p)
            invisible(p)
        })
        names(list) <- genes
        invisible(list)
    }
)



#' @rdname plotMarker
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
        .plotMarkerDimRed(
            object = object,
            genes = genes,
            dimRed = "tsne",
            expression = expression,
            color = color,
            pointsAsNumbers = pointsAsNumbers,
            pointSize = pointSize,
            pointAlpha = pointAlpha,
            label = label,
            labelSize = labelSize,
            dark = dark,
            grid = grid,
            legend = legend,
            aspectRatio = aspectRatio,
            title = title
        )
    }
)



#' @rdname plotMarker
#' @export
setMethod(
    "plotMarkerUMAP",
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
        .plotMarkerDimRed(
            object = object,
            genes = genes,
            dimRed = "umap",
            expression = expression,
            color = color,
            pointsAsNumbers = pointsAsNumbers,
            pointSize = pointSize,
            pointAlpha = pointAlpha,
            label = label,
            labelSize = labelSize,
            dark = dark,
            grid = grid,
            legend = legend,
            aspectRatio = aspectRatio,
            title = title
        )
    }
)



#' @rdname plotMarker
#' @export
setMethod(
    "plotTopMarkers",
    signature("seurat"),
    function(
        object,
        markers,
        dimRed = c("tsne", "umap"),
        headerLevel = 2L,
        ...
    ) {
        validObject(object)
        stopifnot(is(markers, "grouped_df"))
        stopifnot(.isSanitizedMarkers(markers))
        dimRed <- match.arg(dimRed)
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
            plotMarker(
                object = object,
                genes = genes,
                dimRed = dimRed,
                headerLevel = subheaderLevel,
                ...
            )
        })

        invisible(list)
    }
)



#' @rdname plotMarker
#' @export
setMethod(
    "plotKnownMarkersDetected",
    signature("seurat"),
    function(
        object,
        markers,
        dimRed = c("tsne", "umap"),
        headerLevel = 2L,
        ...
    ) {
        assert_has_rows(markers)
        stopifnot(is(markers, "grouped_df"))
        assert_has_rows(markers)
        assert_is_subset("cellType", colnames(markers))
        dimRed <- match.arg(dimRed)
        assertIsAHeaderLevel(headerLevel)

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

            plotMarker(
                object = object,
                genes = genes,
                dimRed = dimRed,
                headerLevel = subheaderLevel,
                ...
            )
        })

        invisible(list)
    }
)
