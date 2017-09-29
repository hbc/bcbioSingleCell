#' Plot Genes per Cell
#'
#' @rdname plotGenesPerCell
#' @name plotGenesPerCell
#' @family Quality Control Metrics
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @inheritParams AllGenerics
#' @param interestingGroup Interesting group, to use for colors.
#' @param min Recommended minimum value cutoff.
#' @param max Recommended maximum value cutoff.
#' @param filterCells Show only the cells that have passed filtering cutoffs.
#'
#' @return [ggplot] grid.
NULL



# Constructors ====
.plotGenesPerCellBoxplot <- function(
    object,
    interestingGroup = "sampleName",
    min = NULL,
    max = NULL,
    filterCells = FALSE) {
    metrics <- metrics(object, filterCells = filterCells)
    medianGenes <- aggregate(nGene ~ sampleID, metrics, median) %>%
        left_join(sampleMetadata(object), by = "sampleID") %>%
        mutate(nGene = round(.data[["nGene"]]))
    p <- ggplot(metrics,
                mapping = aes_string(
                    x = "sampleName",
                    y = "nGene",
                    fill = interestingGroup)) +
        labs(x = "sample",
             y = "genes per cell") +
        geom_boxplot(color = lineColor) +
        geom_label(
            data = medianGenes,
            mapping = aes_(label = ~nGene),
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



.plotGenesPerCellHistogram <- function(
    object,
    min = NULL,
    max = NULL,
    filterCells = FALSE) {
    metrics <- metrics(object, filterCells = filterCells)
    p <- ggplot(metrics,
                aes_(x = ~nGene,
                     fill = ~sampleName)) +
        labs(x = "genes per cell") +
        geom_histogram(bins = bins) +
        scale_x_sqrt() +
        scale_y_sqrt() +
        scale_fill_viridis(discrete = TRUE) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
    if (!is.null(min)) {
        p <- p +
            .qcCutoffLine(xintercept = min)
    }
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



.plotGenesPerCell <- function(
    object,
    interestingGroup,
    min,
    max,
    filterCells = FALSE) {
    if (missing(interestingGroup)) {
        interestingGroup <- interestingGroups(object)[[1]]
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
            max = max,
            filterCells = filterCells),
        .plotGenesPerCellBoxplot(
            object,
            interestingGroup = interestingGroup,
            min = min,
            max = max,
            filterCells = filterCells),
        labels = "auto",
        nrow = 2)
}



# Methods ====
#' @rdname plotGenesPerCell
#' @export
setMethod("plotGenesPerCell", "bcbioSingleCellANY", .plotGenesPerCell)
