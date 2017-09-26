#' Plot Genes per Cell
#'
#' @rdname plotGenesPerCell
#' @name plotGenesPerCell
#' @family Quality Control Metrics
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @param interestingGroup Interesting group, to use for colors.
#' @param min Recommended minimum value cutoff.
#' @param max Recommended maximum value cutoff.
#'
#' @return [ggplot] grid.
NULL



# Constructors ====
.plotGenesPerCellBoxplot <- function(object, interestingGroup, min, max) {
    metrics <- metrics(object)
    medianGenes <- aggregate(nGene ~ sampleID, metrics, median) %>%
        left_join(sampleMetadata(object), by = "sampleID") %>%
        mutate(nGene = round(.data[["nGene"]]))
    p <- ggplot(metrics,
                aes_(x = ~sampleName,
                     y = ~nGene,
                     fill = as.name(interestingGroup))) +
        labs(x = "sample",
             y = "genes per cell") +
        geom_boxplot(color = qcBorderColor) +
        geom_label(
            data = medianGenes,
            aes_(label = ~nGene),
            alpha = qcLabelAlpha,
            color = qcLabelColor,
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
                color = qcCutoffColor,
                linetype = qcLineType,
                size = qcLineSize,
                yintercept = min)
    }
    if (!is.null(max)) {
        p <- p +
            geom_hline(
                alpha = qcLineAlpha,
                color = qcCutoffColor,
                linetype = qcLineType,
                size = qcLineSize,
                yintercept = max)
    }
    if (isTRUE(metadata(object)[["multiplexedFASTQ"]])) {
        p <- p +
            facet_wrap(~fileName)
    }
    p
}



.plotGenesPerCellHistogram <- function(object, min, max) {
    metrics <- metrics(object)
    p <- ggplot(metrics,
                aes_(x = ~nGene,
                     fill = ~sampleName)) +
        labs(x = "genes per cell") +
        geom_histogram(bins = bins) +
        scale_x_sqrt() +
        scale_y_sqrt() +
        scale_fill_viridis(discrete = TRUE) +
        theme(axis.text.x = element_text(angle = 90L, hjust = 1L))
    if (!is.null(min)) {
        p <- p +
            geom_vline(
                alpha = qcLineAlpha,
                color = qcCutoffColor,
                linetype = qcLineType,
                size = qcLineSize,
                xintercept = min)
    }
    if (!is.null(max)) {
        p <- p +
            geom_vline(
                alpha = qcLineAlpha,
                color = qcCutoffColor,
                linetype = qcLineType,
                size = qcLineSize,
                xintercept = max)
    }
    if (isTRUE(metadata(object)[["multiplexedFASTQ"]])) {
        p <- p +
            facet_wrap(~fileName)
    }
    p
}



.plotGenesPerCell <- function(object, interestingGroup, min, max) {
    if (missing(interestingGroup)) {
        interestingGroup <- interestingGroups(object)[[1L]]
    }
    if (missing(min)) {
        min <- object %>%
            metadata() %>%
            .[["filterParams"]] %>%
            .[["minGenes"]]
    }
    if (missing(max)) {
        max <- object %>%
            metadata() %>%
            .[["filterParams"]] %>%
            .[["maxGenes"]]
    }
    plot_grid(
        .plotGenesPerCellHistogram(
            object,
            min = min,
            max = max),
        .plotGenesPerCellBoxplot(
            object,
            interestingGroup = interestingGroup,
            min = min,
            max = max),
        labels = "auto",
        nrow = 2L)
}



# Methods ====
#' @rdname plotGenesPerCell
#' @export
setMethod("plotGenesPerCell", "bcbioSingleCellANY", .plotGenesPerCell)
