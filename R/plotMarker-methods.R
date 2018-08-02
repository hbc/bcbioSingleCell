#' Plot Cell-Type-Specific Gene Markers
#'
#' @description
#' Plot gene expression per cell in multiple formats:
#'
#' 1. [plotMarkerTSNE()]: t-SNE gene expression plot.
#' 2. [plotDot()]: Dot plot.
#' 3. [plotViolin()]: Violin plot.
#'
#' @section Plot top markers:
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
#' @param markers `grouped_df`. Marker genes data.
#'   - [plotTopMarkers()]: must be grouped by "`cluster`".
#'   - [plotKnownMarkersDetected()]: must be grouped by "`cellType`".
#'
#' @return Show graphical output. Invisibly return `ggplot` `list`.
#'
#' @examples
#' # SingleCellExperiment ====
#' object <- cellranger_small
#'
#' # MHC class II genes
#' title <- "MHC class II"
#' genes <- rownames(object)[which(grepl(
#'     pattern = "^major histocompatibility complex, class II",
#'     x = rowData(object)$description
#' ))]
#' print(genes)
#'
#' # t-SNE
#' plotMarkerTSNE(object, genes)
#' plotMarkerTSNE(
#'     object = object,
#'     genes = genes,
#'     expression = "sum",
#'     pointsAsNumbers = TRUE,
#'     dark = TRUE,
#'     label = FALSE,
#'     title = title
#' )
#'
#' # UMAP
#' plotMarkerUMAP(object, genes)
#' plotMarkerTSNE(
#'     object = object,
#'     genes = genes,
#'     expression = "mean",
#'     pointsAsNumbers = TRUE,
#'     dark = TRUE,
#'     label = FALSE,
#'     title = title
#' )
#'
#' # seurat ====
#' object <- seurat_small
#'
#' # Top markers
#' markers <- topMarkers(all_markers_small, n = 1)
#' glimpse(markers)
#' plotTopMarkers(object, markers = tail(markers, 1))
#'
#' # Known markers detected
#' markers <- head(known_markers_small, n = 1)
#' glimpse(markers)
#' plotKnownMarkersDetected(object, markers = head(markers, 1))
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
#' @export
setMethod(
    "plotMarker",
    signature("SingleCellExperiment"),
    function(
        object,
        genes,
        reducedDim = c("TSNE", "UMAP"),
        expression = c("mean", "median", "sum"),
        color = getOption("bcbio.discrete.color", NULL),
        pointSize = getOption("bcbio.pointSize", 0.75),
        pointAlpha = getOption("bcbio.pointAlpha", 0.8),
        pointsAsNumbers = FALSE,
        label = getOption("bcbio.label", TRUE),
        labelSize = getOption("bcbio.labelSize", 6L),
        dark = getOption("bcbio.dark", FALSE),
        grid = getOption("bcbio.grid", FALSE),
        legend = getOption("bcbio.legend", TRUE),
        aspectRatio = getOption("bcbio.aspectRatio", 1L),
        title = TRUE
    ) {
        assert_is_character(genes)
        assert_has_no_duplicates(genes)
        assert_is_subset(genes, rownames(object))
        reducedDim <- match.arg(reducedDim)
        expression <- match.arg(expression)
        # Legacy support for `color = "auto"`
        if (identical(color, "auto")) {
            color <- NULL
        }
        assertIsColorScaleContinuousOrNULL(color)
        assert_is_a_number(pointSize)
        assert_is_a_number(pointAlpha)
        assert_is_a_bool(pointsAsNumbers)
        assert_is_a_bool(label)
        assert_is_a_number(labelSize)
        assert_is_a_bool(dark)
        assert_is_a_bool(legend)
        assert_is_any_of(title, c("character", "logical", "NULL"))
        if (is.character(title)) {
            assert_is_a_string(title)
        }

        # Fetch reduced dimension data
        data <- fetchReducedDimExpressionData(
            object = object,
            genes = genes,
            reducedDim = reducedDim
        )
        axes <- colnames(data)[seq_len(2L)]

        if (!isTRUE(.useGene2symbol(object))) {
            g2s <- gene2symbol(object)
            if (length(g2s)) {
                g2s <- g2s[genes, , drop = FALSE]
                genes <- make.unique(g2s[["geneName"]])
            }
        }
        genes <- sort(unique(genes))

        requiredCols <- c(
            axes,
            "x",
            "y",
            "centerX",
            "centerY",
            "mean",
            "median",
            "ident",
            "sum"
        )
        assert_is_subset(requiredCols, colnames(data))

        p <- ggplot(
            data = data,
            mapping = aes(
                x = !!sym("x"),
                y = !!sym("y"),
                color = !!sym(expression)
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
        p <- p +
            labs(
                x = axes[[1L]],
                y = axes[[2L]],
                title = title,
                subtitle = subtitle
            )

        # Customize legend
        if (isTRUE(legend)) {
            if (is_a_string(genes)) {
                guideTitle <- "expression"
            } else {
                guideTitle <- expression
            }
            p <- p + guides(color = guide_colorbar(title = guideTitle))
        } else {
            p <- p + guides(color = "none")
        }

        if (isTRUE(pointsAsNumbers)) {
            if (pointSize < 4L) pointSize <- 4L
            p <- p +
                geom_text(
                    mapping = aes(
                        x = !!sym("x"),
                        y = !!sym("y"),
                        label = !!sym("ident"),
                        color = !!sym(expression)
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
                    mapping = aes(
                        x = !!sym("centerX"),
                        y = !!sym("centerY"),
                        label = !!sym("ident")
                    ),
                    color = labelColor,
                    size = labelSize,
                    fontface = "bold"
                )
        }

        # Color palette
        if (isTRUE(dark)) {
            theme <- theme_midnight
            if (is.null(color)) {
                color <- darkMarkerColors
            }
        } else {
            theme <- theme_paperwhite
            if (is.null(color)) {
                color <- lightMarkerColors
            }
        }
        p <- p +
            theme(
                aspect_ratio = aspectRatio,
                grid = grid
            )

        if (is(color, "ScaleContinuous")) {
            p <- p + color
        }

        p
    }
)



#' @rdname plotMarker
#' @export
setMethod(
    "plotTopMarkers",
    signature("SingleCellExperiment"),
    function(
        object,
        markers,
        n = 10L,
        direction = c("positive", "negative", "both"),
        coding = FALSE,
        reduction = c("TSNE", "UMAP"),
        headerLevel = 2L,
        ...
    ) {
        # Passthrough: n, direction, coding
        validObject(object)
        stopifnot(is(markers, "grouped_df"))
        stopifnot(.isSanitizedMarkers(markers))
        markers <- topMarkers(
            data = markers,
            n = n,
            direction = direction,
            coding = coding
        )
        reduction <- match.arg(reduction)
        assertIsAHeaderLevel(headerLevel)

        assert_is_subset("cluster", colnames(markers))
        clusters <- levels(markers[["cluster"]])

        list <- pblapply(clusters, function(cluster) {
            genes <- markers %>%
                filter(cluster == !!cluster) %>%
                pull("rowname")
            if (!length(genes)) {
                return(invisible())
            }
            if (length(genes) > 10L) {
                warning("Maximum of 10 genes per cluster is recommended")
            }

            markdownHeader(
                text = paste("Cluster", cluster),
                level = headerLevel,
                tabset = TRUE,
                asis = TRUE
            )

            lapply(genes, function(gene) {
                markdownHeader(
                    text = gene,
                    level = headerLevel + 1L,
                    asis = TRUE
                )
                p <- .plotMarkerReduction(
                    object = object,
                    genes = gene,
                    reduction = reduction,
                    ...
                )
                show(p)
                invisible(p)
            })
        })

        invisible(list)
    }
)



#' @rdname plotMarker
#' @export
setMethod(
    "plotKnownMarkersDetected",
    signature("SingleCellExperiment"),
    function(
        object,
        markers,
        reduction = c("TSNE", "UMAP"),
        headerLevel = 2L,
        ...
    ) {
        assert_has_rows(markers)
        stopifnot(is(markers, "grouped_df"))
        assert_has_rows(markers)
        assert_is_subset("cellType", colnames(markers))
        reduction <- match.arg(reduction)
        assertIsAHeaderLevel(headerLevel)

        cellTypes <- markers %>%
            pull("cellType") %>%
            as.character() %>%
            na.omit() %>%
            unique()
        assert_is_non_empty(cellTypes)

        list <- pblapply(cellTypes, function(cellType) {
            genes <- markers %>%
                filter(cellType == !!cellType) %>%
                pull("geneName") %>%
                as.character() %>%
                na.omit() %>%
                unique()
            assert_is_non_empty(genes)

            markdownHeader(
                text = cellType,
                level = headerLevel,
                tabset = TRUE,
                asis = TRUE
            )

            lapply(genes, function(gene) {
                markdownHeader(
                    text = gene,
                    level = headerLevel + 1L,
                    asis = TRUE
                )
                p <- .plotMarkerReduction(
                    object = object,
                    genes = gene,
                    reduction = reduction,
                    ...
                )
                show(p)
                invisible(p)
            })
        })

        invisible(list)
    }
)
