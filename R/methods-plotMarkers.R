#' Plot Cell-Type Gene Markers
#'
#' @description
#' Plot gene expression per cell in multiple formats:
#'
#' 1. t-SNE plot
#' 2. Violin plot
#' 3. Ridgeline plot
#' 4. Dot plot
#'
#' @name plotMarkers
#' @family Clustering Utilities
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#' @inheritParams plotMarkerTSNE
#' @param tsneColor Color palette to use for tSNE plot.
#' @param violinFill Color palette to use for violin plot.
#' @param dotColor Color palette to use for dot plot.
#'
#' @return No value, only graphical output.
#'
#' @examples
#' load(system.file("extdata/seurat.rda", package = "bcbioSingleCell"))
#' load(system.file("extdata/topMarkers.rda", package = "bcbioSingleCell"))
#'
#' # seurat
#' genes <- head(pull(topMarkers, "geneName"), 2L)
#' plotMarkers(seurat, genes = genes)
NULL



# Constructors =================================================================
#' Plot Marker Seurat Constructor
#'
#' @keywords internal
#' @noRd
#'
#' @param returnAsList Return the `gg` objects as a list.
.plotMarker.seurat <- function(  # nolint
    object,
    gene,
    # FIXME Change tsneColor method here to use function name?
    tsneColor = scale_color_viridis(),
    violinFill = scale_fill_viridis(discrete = TRUE),
    dotColor = scale_color_gradient(
        low = "lightgray",
        high = "purple"
    ),
    dark = TRUE,
    pointsAsNumbers = FALSE,
    title = NULL,
    return = "grid"
) {
    # Parameter integrity ======================================================
    if (!is_a_string(gene)) {
        abort("`gene` must be a string")
    }
    # return
    validReturn <- c("grid", "list")
    if (!return %in% validReturn) {
        abort(paste("`return` must contain:", toString(validReturn)))
    }

    # Plots ====================================================================
    tsne <- plotMarkerTSNE(
        object,
        genes = gene,
        expression = "sum",
        color = tsneColor,
        dark = dark,
        pointsAsNumbers = pointsAsNumbers,
        title = gene,
        subtitle = FALSE
    )
    violin <- plotViolin(
        object,
        genes = gene,
        fill = violinFill
    )
    dot <- plotDot(
        object,
        genes = gene,
        color = dotColor
    )

    # Return ===================================================================
    if (return == "grid") {
        dot <- dot +
            labs(x = "") +
            coord_flip() +
            theme(
                axis.text.y = element_text(angle = 90L, hjust = 0.5),
                legend.position = "none"
            )
        plot_grid(
            tsne,
            violin,
            dot,
            labels = NULL,
            ncol = 1L,
            nrow = 3L,
            rel_heights = c(1L, 0.3, 0.15)
        )
    } else if (return == "list") {
        list(
            "tsne" = tsne,
            "dot" = dot,
            "violin" = violin
        )
    }
}



.plotMarkers.seurat <- function(  # nolint
    object,
    genes,
    # FIXME Change tsneColor method?
    tsneColor = scale_color_viridis(),
    violinFill = scale_fill_viridis(discrete = TRUE),
    dotColor = scale_color_gradient(
        low = "lightgray",
        high = "purple"
    ),
    dark = TRUE,
    pointsAsNumbers = FALSE,
    headerLevel = 2L,
    title = NULL
) {
    list <- lapply(seq_along(genes), function(a) {
        gene <- genes[[a]]
        # Skip and warn if gene is missing
        if (!gene %in% rownames(slot(object, "data"))) {
            return(warn(paste(gene, "missing")))
        }
        if (!is.null(headerLevel)) {
            markdownHeader(gene, level = headerLevel, asis = TRUE)
        }
        p <- .plotMarker.seurat(
            object,
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
        p
    })
    invisible(list)
}



# Methods ======================================================================
#' @rdname plotMarkers
#' @export
setMethod(
    "plotMarkers",
    signature("seurat"),
    .plotMarkers.seurat
)
