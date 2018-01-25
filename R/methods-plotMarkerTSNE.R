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
#' @param format Gene identifier format. Defaults to `symbol` but `ensgene`
#'   is also supported, as long as Ensembl gene to symbol identifier mappings
#'   are defined.
#' @param colorPoints Color points by geometric mean (`geomean`) or expression
#'   of individual gene (`expression`).
#' @param legend Show plot legend.
#' @param subtitle Include gene(s) in the subtitle.
#'
#' @return [ggplot].
#'
#' @examples
#' load(system.file(
#'     file.path("extdata", "seurat.rda"),
#'     package = "bcbioSingleCell"))
#'
#' symbol <- counts(seurat) %>% rownames() %>% .[[1]]
#' print(symbol)
#'
#' ensgene <- bcbio(seurat, "gene2symbol") %>%
#'     .[which(.[["symbol"]] %in% symbol), "ensgene", drop = TRUE]
#' print(ensgene)
#'
#' # seurat
#' plotMarkerTSNE(seurat, genes = symbol, format = "symbol")
#' plotMarkerTSNE(seurat, genes = ensgene, format = "ensgene")
#'
#' # data.frame
#' df <- fetchTSNEExpressionData(seurat, genes = symbol)
#' plotMarkerTSNE(df)
NULL



# Constructors =================================================================
#' Plot Marker tSNE Constructor
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom bcbioBase midnightTheme
#' @importFrom ggplot2 aes_string geom_point geom_text ggplot labs guides
#'   guide_colorbar theme
#' @importFrom viridis scale_color_viridis
#'
#' @param object Marker gene expression [data.frame] returned from
#'   [fetchTSNEExpressionData()].
.plotMarkerTSNE <- function(
    object,
    colorPoints = "geomean",
    pointsAsNumbers = FALSE,
    pointSize = 0.5,
    pointAlpha = 0.8,
    label = TRUE,
    labelSize = 6L,
    color = viridis::scale_color_viridis(),
    dark = TRUE,
    legend = TRUE,
    title = NULL,
    subtitle = TRUE) {
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
        abort(paste(
            "Required columns:", toString(requiredCols)
        ))
    }
    validColorPoints <- c("expression", "geomean")
    if (!colorPoints %in% validColorPoints) {
        abort(paste(
            "`colorPoints` must contain:",
            toString(validColorPoints)
        ))
    }

    # Prepare a list of the genes used for the ggplot subtitle
    genes <- unique(pull(object, "gene"))
    # Use `expression` if we're only plotting a single gene. The `geomean`
    # argument for `colorPoints` is only informative for 2+ genes.
    if (length(genes) == 1L) {
        colorPoints <- "expression"
    }
    # Limit to displaying the top 5, with an ellipsis if necessary
    if (length(genes) > 5L) {
        genes <- c(genes[1L:5L], "...")
    }
    genes <- toString(genes)

    p <- ggplot(
        object,
        mapping = aes_string(
            x = "tSNE1",
            y = "tSNE2",
            color = colorPoints)
    )

    if (isTRUE(dark)) {
        p <- p + midnightTheme()
    }

    # Labels
    if (!is.character(title)) {
        title <- NULL
    }
    if (isTRUE(subtitle)) {
        subtitle <- genes
    } else {
        subtitle <- NULL
    }
    p <- p +
        labs(title = title,
             subtitle = subtitle)

    # Customize legend
    if (isTRUE(legend)) {
        p <- p +
            # Make the guide longer than normal, to improve appearance of values
            # containing a decimal point
            guides(color = guide_colorbar(
                barwidth = 20L,
                barheight = 1L,
                direction = "horizontal")) +
            theme(legend.justification = "center",
                  legend.position = "bottom")
    } else {
        p <- p + theme(legend.position = "none")
    }

    if (isTRUE(pointsAsNumbers)) {
        p <- p +
            geom_text(
                mapping = aes_string(
                    x = "tSNE1",
                    y = "tSNE2",
                    label = "ident",
                    color = colorPoints),
                alpha = pointAlpha,
                size = pointSize)
    } else {
        p <- p +
            geom_point(
                alpha = pointAlpha,
                size = pointSize)
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
                size = labelSize,
                fontface = "bold")
    }

    if (is(color, "ScaleContinuous")) {
        p <- p + color
    }

    p
}



# Methods ======================================================================
#' @rdname plotMarkerTSNE
#' @export
setMethod(
    "plotMarkerTSNE",
    signature("data.frame"),
    .plotMarkerTSNE)



#' @rdname plotMarkerTSNE
#' @export
setMethod(
    "plotMarkerTSNE",
    signature("seurat"),
    function(
        object,
        genes,
        format = "symbol",
        colorPoints = "geomean",
        pointsAsNumbers = FALSE,
        pointSize = 0.5,
        pointAlpha = 0.8,
        label = TRUE,
        labelSize = 6L,
        color = viridis::scale_color_viridis(),
        dark = TRUE,
        title = NULL,
        subtitle = TRUE) {
        .checkFormat(format)
        if (format == "ensgene") {
            genes <- .convertGenesToSymbols(object, genes = genes)
        }
        data <- fetchTSNEExpressionData(object, genes = genes)
        .plotMarkerTSNE(
            object = data,
            colorPoints = colorPoints,
            pointsAsNumbers = pointsAsNumbers,
            pointSize = pointSize,
            pointAlpha = pointAlpha,
            label = label,
            labelSize = labelSize,
            color = color,
            dark = dark,
            title = title,
            subtitle = subtitle)
    })
