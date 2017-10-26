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
#'
#' @param color Color palette to use for points. Defaults to the inferno
#'   palette from the viridis package.
#' @param colorPoints Color points by geometric mean (`geomean`) or expression
#'   of individual gene (`expression`).
#' @param dark Enable dark mode.
#' @param pointsAsNumbers Plot the points as numbers (`TRUE`) or dots (`FALSE`).
#' @param title Plot title.
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
#'     dplyr::filter(cluster == 0) %>%
#'     dplyr::pull(symbol)
#'
#' # Fetch the t-SNE expression data for the desired gene symbols
#' dat <- fetchTSNEExpressionData(seurat, genes = genes)
#' print(unique(dat$gene))
#'
#' # To make t-SNE plot colored by geometric mean of topGenes
#' plotTSNE(dat, colorPoints = "geomean")
#'
#' # To make faceted t-SNE plot of each gene (looks good at up to 6 genes)
#' plotTSNE(dat, colorPoints = "expression") +
#'     ggplot2::facet_wrap(~gene)
#' }
NULL



# Constructors ====
#' Plot Marker tSNE Constructor
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom ggplot2 aes_string geom_point geom_text ggplot labs guides
#'   guide_colorbar
#' @importFrom viridis scale_color_viridis
#'
#' @param object Marker gene expression [tibble].
.plotMarkerTSNE <- function(
    object,
    colorPoints = "geomean",
    color = scale_color_viridis(option = "inferno"),
    dark = TRUE,
    pointsAsNumbers = FALSE,
    label = TRUE,
    title = NULL) {
    if (!colorPoints %in% c("expression", "geomean")) {
        stop("colorPoints supports 'geomean' or 'expression'", call. = FALSE)
    }
    # Check to make sure only a subset of markers is passed in
    genes <- sort(unique(object[["gene"]]))
    if (length(genes) > 100) {
        stop("geomean argument supports up to 100 marker genes", call. = FALSE)
    }
    p <- object %>%
        ggplot(
            mapping = aes_string(
                x = "tSNE1",
                y = "tSNE2",
                color = colorPoints)
        ) +
        labs(title = title,
             subtitle = genes) +
        guides(color = guide_colorbar(direction = "horizontal"))
    if (isTRUE(pointsAsNumbers)) {
        # This seems to take longer to plot than the points?
        p <- p +
            geom_text(
                mapping = aes_string(
                    x = "tSNE1",
                    y = "tSNE2",
                    label = "ident",
                    color = colorPoints),
                size = 2)
    } else {
        p <- p + geom_point()
    }
    if (isTRUE(label)) {
        p <- p +
            geom_text(
                mapping = aes_string(
                    x = "centerX",
                    y = "centerY",
                    label = "ident"),
                color = "white",
                size = 6,
                fontface = "bold")
    }
    if (!is.null(color)) {
        p <- p + color
    }
    if (isTRUE(dark)) {
        p <- p + darkTheme()
    }
    p
}



# Methods ====
#' @rdname plotMarkerTSNE
#' @export
setMethod("plotMarkerTSNE", "grouped_df", function(
    object,
    colorPoints = "geomean",
    color = scale_color_viridis(option = "inferno"),
    dark = TRUE,
    pointsAsNumbers = FALSE,
    label = TRUE) {
    requiredCols <- c(
        "centerX",
        "centerY",
        "expression",
        "gene",
        "geomean",
        "ident",
        "tSNE1",
        "tSNE2")
    if (!all(requiredCols %in% colnames(object))) {
        stop(paste(
            "Required columns:", toString(requiredCols)
        ), call. = FALSE)
    }
    .plotMarkerTSNE(
        object = object,
        colorPoints = colorPoints,
        color = color,
        dark = dark,
        pointsAsNumbers = pointsAsNumbers,
        label = label)
})



#' @rdname plotMarkerTSNE
#' @export
setMethod("plotMarkerTSNE", "seurat", function(
    object,
    genes,
    colorPoints = "geomean",
    color = scale_color_viridis(option = "inferno"),
    dark = TRUE,
    pointsAsNumbers = FALSE,
    label = TRUE,
    title = NULL) {
    fetch <- fetchTSNEExpressionData(object, genes = genes)
        .plotMarkerTSNE(
            object = fetch,
            colorPoints = colorPoints,
            color = color,
            dark = dark,
            pointsAsNumbers = pointsAsNumbers,
            label = label,
            title = title)
})
