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
.plotNoveltyBoxplot <- function(object, interestingGroup, min) {
    metrics <- metrics(object)
    medianNovelty <-
        aggregate(log10GenesPerUMI ~ sampleID, metrics, median) %>%
        left_join(sampleMetadata(object), by = "sampleID")
    p <- ggplot(metrics,
                aes_(x = ~sampleName,
                     y = ~log10GenesPerUMI,
                     fill = as.name(interestingGroup))) +
        labs(x = "sample",
             y = "log10 genes per UMI") +
        geom_boxplot(colour = lineColor) +
        geom_label(
            data = medianNovelty,
            aes_(label = ~round(log10GenesPerUMI, digits = 2L)),
            alpha = qcLabelAlpha,
            colour = qcLabelColor,
            fill = qcLabelFill,
            fontface = qcLabelFontface,
            label.padding = qcLabelPadding,
            label.size = qcLabelSize,
            show.legend = FALSE) +
        scale_y_sqrt() +
        scale_fill_viridis(discrete = TRUE) +
        theme(axis.text.x = element_text(angle = 90L, hjust = 1L))
    if (!is.null(min)) {
        p <- p +
            geom_hline(
                alpha = qcLineAlpha,
                colour = qcCutoffColor,
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



.plotNoveltyHistogram <- function(object, min) {
    metrics <- metrics(object)
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
                colour = qcCutoffColor,
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



.plotNovelty <- function(object, interestingGroup, min) {
    if (missing(interestingGroup)) {
        interestingGroup <- interestingGroups(object)[[1L]]
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
            min = min),
        .plotNoveltyBoxplot(
            object,
            interestingGroup = interestingGroup,
            min = min),
        labels = "auto",
        nrow = 2L)
}



# Methods ====
#' @rdname plotNovelty
#' @export
setMethod("plotNovelty", "bcbioSingleCellANY", .plotNovelty)
