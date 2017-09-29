#' Plot Mitochondrial Transcript Abundance
#'
#' @rdname plotMitoRatio
#' @name plotMitoRatio
#' @family Quality Control Metrics
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @inherit plotGenesPerCell
NULL



# Constructors ====
.plotMitoRatioBoxplot <- function(
    object,
    interestingGroup = "sampleName",
    max = NULL,
    filterCells = FALSE) {
    metrics <- metrics(object, filterCells = filterCells)
    medianMitoRatio <-
        aggregate(mitoRatio ~ sampleID, metrics, median) %>%
        left_join(sampleMetadata(object), by = "sampleID")
    p <- ggplot(metrics,
                aes_(x = ~sampleName,
                     y = ~mitoRatio,
                     fill = as.name(interestingGroup))) +
        geom_boxplot(color = lineColor) +
        geom_label(
            data = medianMitoRatio,
            aes_(label = ~round(mitoRatio, digits = 2)),
            alpha = qcLabelAlpha,
            color = qcLabelColor,
            fill = qcLabelFill,
            fontface = qcLabelFontface,
            label.padding = qcLabelPadding,
            label.size = qcLabelSize,
            show.legend = FALSE) +
        scale_y_sqrt() +
        scale_fill_viridis(discrete = TRUE) +
        labs(x = "sample",
             y = "mito ratio") +
        theme(axis.text.x = element_text(angle = 15, hjust = 1))
    if (!is.null(max)) {
        p <- p +
            .qcCutoffLine(yintercept = max)
    }
    if (isTRUE(metadata(object)[["multiplexedFASTQ"]])) {
        p <- p +
            facet_wrap(~fileName)
    }
    p
}



.plotMitoRatioHistogram <- function(
    object,
    max = NULL,
    filterCells = FALSE) {
    metrics <- metrics(object, filterCells = filterCells)
    p <- ggplot(metrics,
                aes_(x = ~mitoRatio,
                     fill = ~sampleName)) +
        geom_histogram(bins = bins) +
        scale_x_sqrt() +
        scale_y_sqrt() +
        scale_fill_viridis(discrete = TRUE) +
        labs(x = "mito ratio",
             y = "sample")
    if (!is.null(max)) {
        p <- p +
            .qcCutoffLine(xintercept = max)
    }
    if (isTRUE(metadata(object)[["multiplexedFASTQ"]])) {
        p <- p +
            facet_wrap(~fileName)
    }
    p
}



.plotMitoRatioScatterplot <- function(
    object,
    filterCells = FALSE) {
    metrics <- metrics(object, filterCells = filterCells)
    p <- ggplot(metrics,
                aes_(x = ~nCoding,
                     y = ~nMito,
                     color = ~sampleName)) +
        labs(x = "mito counts",
             y = "coding counts") +
        geom_point(alpha = 0.6, size = 1) +
        geom_smooth(se = FALSE) +
        scale_x_sqrt() +
        scale_y_sqrt() +
        scale_color_viridis(discrete = TRUE)
    if (isTRUE(metadata(object)[["multiplexedFASTQ"]])) {
        p <- p +
            facet_wrap(~fileName)
    }
    p
}



.plotMitoRatio <- function(
    object,
    interestingGroup,
    max,
    filterCells = FALSE) {
    if (missing(interestingGroup)) {
        interestingGroup <- interestingGroups(object)[[1]]
    }
    if (missing(max)) {
        max <- object %>%
            metadata() %>%
            .[["filterParams"]] %>%
            .[["maxMitoRatio"]]
    }
    plot_grid(
        .plotMitoRatioScatterplot(object, filterCells = filterCells),
        .plotMitoRatioHistogram(
            object,
            max = max,
            filterCells = filterCells),
        .plotMitoRatioBoxplot(
            object,
            interestingGroup = interestingGroup,
            max = max,
            filterCells = filterCells),
        labels = "auto",
        nrow = 3)
}



# Methods ====
#' @rdname plotMitoRatio
#' @export
setMethod("plotMitoRatio", "bcbioSingleCellANY", .plotMitoRatio)
