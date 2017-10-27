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
#'   palette from the viridis package. Use [ggplot2::scale_color_gradient()]
#'   to easily define your own low/high gradient.
#' @param colorPoints Color points by geometric mean (`geomean`) or expression
#'   of individual gene (`expression`).
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
#' @param fetchData Marker gene expression from [fetchTSNEExpressionData()].
.plotMarkerTSNE <- function(
    fetchData,
    colorPoints = "geomean",
    color = scale_color_viridis(option = "inferno"),
    dark = TRUE,
    pointsAsNumbers = FALSE,
    label = TRUE,
    title = NULL) {
    if (!colorPoints %in% c("expression", "geomean")) {
        stop("colorPoints supports 'geomean' or 'expression'", call. = FALSE)
    }
    # Prepare a list of the genes used for the ggplot subtitle
    genes <- pull(fetchData, "gene") %>%
        unique()
    # Use `expression` if we're only plotting a single gene. The `geomean`
    # argument for `colorPoints` is only informative for 2+ genes.
    if (length(genes) == 1) {
        colorPoints = "expression"
    }
    # Limit to displaying the top 5, with an ellipsis if necessary
    if (length(genes) > 5) {
        genes <- c(genes[1:5], "...")
    }
    genes <- toString(genes)
    p <- fetchData %>%
        ggplot(
            mapping = aes_string(
                x = "tSNE1",
                y = "tSNE2",
                color = colorPoints)
        )
    if (isTRUE(dark)) {
        p <- p + darkTheme()
    }
    p <- p +
        labs(title = title,
         subtitle = genes) +
        # Make the guide longer than normal, to improve appearance of values
        # containing a decimal point
        guides(color = guide_colorbar(
            barwidth = 20,
            barheight = 1,
            direction = "horizontal")) +
        theme(legend.justification = "center",
              legend.position = "bottom")
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
                    label = "ident"),
                color = labelColor,
                size = 6,
                fontface = "bold")
    }
    if (!is.null(color)) {
        p <- p + color
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
        fetchData = object,
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
    fetchData <- fetchTSNEExpressionData(object, genes = genes)
        .plotMarkerTSNE(
            fetchData = fetchData,
            colorPoints = colorPoints,
            color = color,
            dark = dark,
            pointsAsNumbers = pointsAsNumbers,
            label = label,
            title = title)
})
