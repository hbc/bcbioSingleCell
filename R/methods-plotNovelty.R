#' Plot Novelty Score
#'
#' @rdname plotNovelty
#' @name plotNovelty
#' @family Quality Control Metrics
#' @author Michael Steinbaugh
#'
#' @inherit plotGenesPerCell
#'
#' @details "Novelty" refers to log10 genes detected per count.
NULL



# Constructors ====
.plotNoveltyBoxplot <- function(
    object,
    interestingGroup = "sampleName",
    min = NULL,
    filterCells = FALSE) {
    metrics <- metrics(object, filterCells = filterCells)
    medianNovelty <-
        aggregate(log10GenesPerUMI ~ sampleID, metrics, median) %>%
        left_join(sampleMetadata(object), by = "sampleID")
    p <- ggplot(metrics,
                aes_(x = ~sampleName,
                     y = ~log10GenesPerUMI,
                     fill = as.name(interestingGroup))) +
        labs(x = "sample",
             y = "log10 genes per UMI") +
        geom_boxplot(color = lineColor) +
        geom_label(
            data = medianNovelty,
            aes_(label = ~round(log10GenesPerUMI, digits = 2)),
            alpha = qcLabelAlpha,
            color = qcLabelColor,
            fill = qcLabelFill,
            fontface = qcLabelFontface,
            label.padding = qcLabelPadding,
            label.size = qcLabelSize,
            show.legend = FALSE) +
        scale_y_sqrt() +
        scale_fill_viridis(discrete = TRUE) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
    if (!is.null(min)) {
        p <- p +
            geom_hline(
                alpha = qcLineAlpha,
                color = qcCutoffColor,
                linetype = qcLineType,
                size = qcLineSize,
                yintercept = min)
    }
    if (isTRUE(metadata(object)[["multiplexedFASTQ"]])) {
        p <- p +
            facet_wrap(~fileName)
    }
    p
}



.plotNoveltyHistogram <- function(
    object,
    min = NULL,
    filterCells = FALSE) {
    if (isTRUE(filterCells)) {
    }
    metrics <- metrics(object, filterCells = filterCells)
    p <- ggplot(metrics,
                aes_(x = ~log10GenesPerUMI,
                     fill = ~sampleName)) +
        labs(x = "log10 genes per UMI") +
        geom_histogram(bins = bins) +
        scale_x_sqrt() +
        scale_y_sqrt() +
        scale_fill_viridis(discrete = TRUE)
    if (!is.null(min)) {
        p <- p +
            geom_vline(
                alpha = qcLineAlpha,
                color = qcCutoffColor,
                linetype = qcLineType,
                size = qcLineSize,
                xintercept = min)
    }
    if (isTRUE(metadata(object)[["multiplexedFASTQ"]])) {
        p <- p +
            facet_wrap(~fileName)
    }
    p
}



.plotNovelty <- function(
    object,
    interestingGroup,
    min,
    filterCells = FALSE) {
    if (missing(interestingGroup)) {
        interestingGroup <- interestingGroups(object)[[1]]
    }
    if (missing(min)) {
        min <- object %>%
            metadata() %>%
            .[["filterParams"]] %>%
            .[["minNovelty"]]
    }
    plot_grid(
        .plotNoveltyHistogram(
            object,
            min = min,
            filterCells = filterCells),
        .plotNoveltyBoxplot(
            object,
            interestingGroup = interestingGroup,
            min = min,
            filterCells = filterCells),
        labels = "auto",
        nrow = 2)
}



# Methods ====
#' @rdname plotNovelty
#' @export
setMethod("plotNovelty", "bcbioSingleCellANY", .plotNovelty)
