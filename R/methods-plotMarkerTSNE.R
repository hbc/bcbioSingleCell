#' Plot Marker t-SNE
#'
#' @rdname plotMarkerTSNE
#' @name plotMarkerTSNE
#' @family t-SNE Utilities
#' @author Rory Kirchner, Michael Steinbaugh
#'
#' @inheritParams AllGenerics
#' @inheritParams fetchTSNEExpressionData
#' @inheritParams plotTSNE
#' @param colorpoints Color points by geometric mean or expression of
#'  individual gene.
#' @param pointsAsNumbers Plot the points as numbers (`TRUE`) or dots (`FALSE`).
#'
#' @return [ggplot].
#'
#' @examples
#' \dontrun{
#' data(markers, seurat)
#' top <- topMarkers(markers)
#'
#' # Let's take the top markers specific to cluster 0, as an example
#' genes <- top %>%
#'     filter(cluster == 0) %>%
#'     pull(symbol)
#'
#' # Fetch the t-SNE expression data for the desired gene symbols
#' dat <- fetchTSNEExpressionData(seurat, genes = genes)
#' print(unique(dat$gene))
#'
#' # To make t-SNE plot colored by geometric mean of topGenes
#' plotTSNE(dat, colorpoints = "geomean")
#'
#' # To make faceted t-SNE plot of each gene (looks good at up to 6 genes)
#' plotTSNE(dat, colorpoints = "expression") +
#'     facet_wrap(~gene)
#' }
NULL



# Constructors ====
.plotMarkerTSNE <- function(
    tibble,
    colorpoints = "geomean",
    pointsAsNumbers = FALSE,
    label = TRUE) {
    if (!colorpoints %in% c("expression", "geomean")) {
        stop("colorpoints supports 'geomean' or 'expression'")
    }
    # Check to make sure only a subset of markers is passed in
    genes <- sort(unique(tibble[["gene"]]))
    if (length(genes) > 100) {
        stop("geomean argument supports up to 100 marker genes")
    }
    p <- tibble %>%
        ggplot(
            aes_string(
                x = "tSNE_1",
                y = "tSNE_2",
                color = colorpoints))
    if (isTRUE(pointsAsNumbers)) {
        # This seems to take longer to plot than the points?
        p <- p +
            geom_text(
                aes_string(
                    x = "tSNE_1",
                    y = "tSNE_2",
                    label = "object.ident",
                    color = colorpoints),
                size = 2)
    } else {
        p <- p +
            geom_point()
    }
    if (isTRUE(label)) {
        p <- p +
            geom_text(
                aes_string(
                    x = "centerx",
                    y = "centery",
                    label = "object.ident"),
                color = "white",
                size = 6,
                fontface = "bold")
    }
    p <- p +
        DarkTheme() +
        scale_color_viridis(option = "inferno") +
        guides(color = guide_colorbar(direction = "horizontal")) +
        theme(legend.justification = "center",
              legend.position = "bottom")
    if (length(genes) <= 10) {
        p <- p +
            ggtitle(toString(genes))
    }
    p
}



# Methods ====
#' @rdname plotMarkerTSNE
#' @export
setMethod("plotMarkerTSNE", "grouped_df", function(
    object,
    colorpoints = "geomean",
    pointsAsNumbers = FALSE,
    label = TRUE) {
    requiredCols <- c(
        "centerx",
        "centery",
        "expression",
        "gene",
        "geomean",
        "object.ident",
        "tSNE_1",
        "tSNE_2")
    if (!all(requiredCols %in% colnames(object))) {
        stop(paste("Required columns:", toString(requiredCols)),
             call. = FALSE)
    }
    .plotMarkerTSNE(
        tibble = object,
        colorpoints = colorpoints,
        pointsAsNumbers = pointsAsNumbers,
        label = label)
})



#' @rdname plotMarkerTSNE
#' @export
setMethod("plotMarkerTSNE", "seurat", function(
    object,
    genes,
    colorpoints = "geomean",
    pointsAsNumbers = FALSE,
    label = TRUE) {
    tibble <- fetchTSNEExpressionData(object, genes = genes)
        .plotMarkerTSNE(
            tibble = tibble,
            colorpoints = colorpoints,
            pointsAsNumbers = pointsAsNumbers,
            label = label)
})
