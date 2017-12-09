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
#' @importFrom basejump midnightTheme
#' @importFrom ggplot2 aes_string geom_point geom_text ggplot labs guides
#'   guide_colorbar theme
#' @importFrom viridis scale_color_viridis
#'
#' @param data Marker gene expression from [fetchTSNEExpressionData()].
.plotMarkerTSNE <- function(
    data,
    colorPoints = "geomean",
    pointsAsNumbers = FALSE,
    pointSize = 1,
    label = TRUE,
    labelSize = 6,
    color = scale_color_viridis(),
    dark = TRUE,
    legend = TRUE,
    title = NULL) {
    if (!colorPoints %in% c("expression", "geomean")) {
        stop("colorPoints supports 'geomean' or 'expression'", call. = FALSE)
    }
    # Prepare a list of the genes used for the ggplot subtitle
    genes <- unique(pull(data, "gene"))
    # Use `expression` if we're only plotting a single gene. The `geomean`
    # argument for `colorPoints` is only informative for 2+ genes.
    if (length(genes) == 1) {
        colorPoints <- "expression"
    }
    # Limit to displaying the top 5, with an ellipsis if necessary
    if (length(genes) > 5) {
        genes <- c(genes[1:5], "...")
    }
    genes <- toString(genes)
    p <- ggplot(
        data,
        mapping = aes_string(
            x = "tSNE1",
            y = "tSNE2",
            color = colorPoints)
    )
    if (isTRUE(dark)) {
        p <- p + midnightTheme()
    }
    if (isTRUE(legend)) {
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
                size = pointSize)
    } else {
        p <- p + geom_point(size = pointSize)
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
    if (!is.null(color)) {
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
    function(
        object,
        colorPoints = "geomean",
        pointsAsNumbers = FALSE,
        pointSize = 1,
        label = TRUE,
        labelSize = 6,
        color = scale_color_viridis(),
        dark = TRUE,
        legend = TRUE,
        title = NULL) {
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
            data = object,
            colorPoints = colorPoints,
            pointsAsNumbers = pointsAsNumbers,
            pointSize = pointSize,
            label = label,
            labelSize = labelSize,
            color = color,
            dark = dark,
            legend = legend,
            title = title)
    })



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
        pointSize = 1,
        label = TRUE,
        labelSize = 6,
        color = scale_color_viridis(),
        dark = TRUE,
        title = NULL) {
        validFormat <- c("ensgene", "symbol")
        if (!format %in% validFormat) {
            stop(paste(
                "'format' must contain:", toString(validFormat)
            ), call. = FALSE)
        }
        rownames <- rownames(slot(object, "data"))
        if (format == "ensgene") {
            ensgene <- genes
            gene2symbol <- bcbio(object, "gene2symbol")
            if (is.null(gene2symbol)) {
                stop(paste(
                    "seurat object doesn't contain",
                    "Ensembl gene to symbol mappings"
                ), call. = FALSE)
            }
            genes <- gene2symbol %>%
                .[which(.[["ensgene"]] %in% ensgene), "symbol", drop = TRUE]
        }
        data <- fetchTSNEExpressionData(object, genes = genes)
        .plotMarkerTSNE(
            data = data,
            colorPoints = colorPoints,
            pointsAsNumbers = pointsAsNumbers,
            pointSize = pointSize,
            label = label,
            labelSize = labelSize,
            color = color,
            dark = dark,
            title = title)
    })
