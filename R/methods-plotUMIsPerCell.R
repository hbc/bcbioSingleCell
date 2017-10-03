#' Plot UMIs per Cell
#'
#' Plot the universal molecular identifiers (UMIs) per cell.
#'
#' @rdname plotUMIsPerCell
#' @name plotUMIsPerCell
#' @family Quality Control Metrics
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @inherit plotGenesPerCell
NULL



# Constructors ====
.plotUMIsPerCellBoxplot <- function(
    object,
    interestingGroup = "sampleName",
    min = NULL,
    filterCells = FALSE) {
    metrics <- metrics(object, filterCells = filterCells)
    p <- ggplot(metrics,
                mapping = aes_string(
                    x = "sampleName",
                    y = "nUMI",
                    fill = interestingGroup)) +
        geom_boxplot(color = lineColor) +
        scale_y_log10() +
        scale_fill_viridis(discrete = TRUE) +
        labs(x = "sample",
             y = "umis per cell") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
    if (length(unique(metrics[["sampleName"]])) <= qcLabelMaxNum) {
        medianUMIs <- aggregate(nUMI ~ sampleID, metrics, median) %>%
            left_join(sampleMetadata(object), by = "sampleID") %>%
            mutate(nUMI = round(.data[["nUMI"]]))
        p <- p +
            geom_label(
                data = medianUMIs,
                mapping = aes_string(label = "nUMI"),
                alpha = qcLabelAlpha,
                color = qcLabelColor,
                fill = qcLabelFill,
                fontface = qcLabelFontface,
                label.padding = qcLabelPadding,
                label.size = qcLabelSize,
                show.legend = FALSE)
    }
    if (!is.null(min)) {
        p <- p +
            .qcCutoffLine(yintercept = min)
    }
    if (isTRUE(metadata(object)[["multiplexedFASTQ"]])) {
        p <- p +
            facet_wrap(~fileName)
    }
    p
}



.plotUMIsPerCellHistogram <- function(
    object,
    min = NULL,
    filterCells = FALSE) {
    metrics <- metrics(object, filterCells = filterCells)
    p <- ggplot(metrics,
                mapping = aes_string(
                    x = "nUMI",
                    fill = "sampleName")) +
        labs(x = "umis per cell") +
        geom_histogram(bins = bins) +
        scale_x_log10() +
        scale_y_sqrt() +
        scale_fill_viridis(discrete = TRUE) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
    if (!is.null(min)) {
        p <- p +
            .qcCutoffLine(xintercept = min)
    }
    if (isTRUE(metadata(object)[["multiplexedFASTQ"]])) {
        p <- p +
            facet_wrap(~fileName)
    }
    p
}



.plotUMIsPerCell <- function(
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
            .[["minUMIs"]]
    }
    plot_grid(
        .plotUMIsPerCellHistogram(
            object,
            min = min,
            filterCells = filterCells),
        .plotUMIsPerCellBoxplot(
            object,
            interestingGroup = interestingGroup,
            min = min,
            filterCells = filterCells),
        labels = "auto",
        nrow = 2
    )
}



# Methods ====
#' @rdname plotUMIsPerCell
#' @export
setMethod("plotUMIsPerCell", "bcbioSingleCell", .plotUMIsPerCell)
