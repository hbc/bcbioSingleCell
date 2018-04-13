#' Plot Cell-Type Gene Marker(s)
#'
#' @description
#' Plot gene expression per cell in multiple formats:
#'
#' 1. t-SNE plot
#' 2. Violin plot
#' 3. Ridgeline plot
#' 4. Dot plot
#'
#' @name plotMarker
#' @family Clustering Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#' @param tsneColor Color palette to use for tSNE plot.
#' @param violinFill Color palette to use for violin plot.
#' @param dotColor Color palette to use for dot plot.
#'
#' @return Show graphical output. Invisibly return `ggplot` `list`.
#'
#' @examples
#' # seurat ===
#' top <- topMarkers(all_markers_small, n = 1L)
#' genes <- pull(top, "rowname")
#'
#' # Individual gene
#' gene <- genes[[1L]]
#' plotMarker(seurat_small, gene = gene, dark = TRUE)
#' plotMarker(seurat_small, gene = gene, dark = FALSE)
#'
#' # Multiple genes
#' plotMarkers(seurat_small, genes = genes)
NULL



# Methods ======================================================================
#' @rdname plotMarker
#' @export
setMethod(
    "plotMarker",
    signature("seurat"),
    function(
        object,
        gene,
        dark = FALSE,
        return = c("grid", "list"),
        ...
    ) {
        assert_is_a_string(gene)
        assert_is_a_bool(dark)
        return <- match.arg(return)

        # Plots ================================================================
        if (isTRUE(dark)) {
            violinFill <- "white"
            dotColor <- scale_color_gradient(
                low = "black",
                high = "white"
            )
        } else {
            violinFill <- "black"
            dotColor <- scale_color_gradient(
                low = "white",
                high = "black"
            )
        }

        tsne <- plotMarkerTSNE(
            object,
            genes = gene,
            expression = "sum",
            dark = dark,
            legend = TRUE,
            title = gene,
            subtitle = FALSE,
            ...
        )

        dot <- plotDot(
            object,
            genes = gene,
            color = dotColor,
            dark = dark
        )

        violin <- plotViolin(
            object,
            genes = gene,
            scale = "width",
            fill = violinFill,
            dark = dark,
            return = "list"
        ) %>%
            # Get the ggplot object from the list return
            .[[1L]]

        # Return ===============================================================
        if (return == "grid") {
            violin <- violin +
                .minimalAxes()
            dot <- dot +
                coord_flip() +
                .minimalAxes()
            p <- plot_grid(
                tsne,
                dot,
                violin,
                labels = NULL,
                ncol = 1L,
                nrow = 3L,
                rel_heights = c(1L, 0.1, 0.1)
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
setMethod(
    "plotMarkers",
    signature("seurat"),
    function(
        object,
        genes,
        tsneColor,
        violinFill = scale_fill_viridis(discrete = TRUE),
        dotColor = scale_color_gradient(
            low = "lightgray",
            high = "orange"
        ),
        dark = FALSE,
        pointsAsNumbers = FALSE,
        headerLevel = 2L,
        title = NULL
    ) {
        assert_is_character(genes)
        assert_is_subset(genes, rownames(object))
        assertIsAHeaderLevel(headerLevel)

        list <- lapply(genes, function(gene) {
            markdownHeader(gene, level = headerLevel, asis = TRUE)
            p <- plotMarker(
                object = object,
                gene = gene,
                tsneColor = tsneColor,
                violinFill = violinFill,
                dotColor = dotColor,
                dark = dark,
                pointsAsNumbers = pointsAsNumbers,
                title = title,
                return = "grid"
            )
            show(p)
            invisible(p)
        })

        invisible(list)
    }
)
