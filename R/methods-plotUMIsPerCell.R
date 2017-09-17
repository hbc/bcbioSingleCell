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
.plotUMIsPerCellBoxplot <- function(object, min) {
    metrics <- metrics(object)
    medianUMIs <- aggregate(nUMI ~ sampleID, metrics, median) %>%
        left_join(sampleMetadata(object), by = "sampleID") %>%
        mutate(nUMI = round(.data[["nUMI"]]))
    interestingGroup <- interestingGroups(object)[[1L]]
    p <- ggplot(metrics,
        aes_(x = ~sampleName,
             y = ~nUMI,
             fill = as.name(interestingGroup))) +
        labs(x = "sample",
             y = "umis per cell") +
        geom_boxplot() +
        geom_label(data = medianUMIs,
                   aes_(label = ~nUMI),
                   alpha = qcLabelAlpha,
                   color = qcLabelColor,
                   fill = qcLabelFill,
                   fontface = qcLabelFontface,
                   label.padding = qcLabelPadding,
                   label.size = qcLabelSize,
                   show.legend = FALSE) +
        geom_hline(alpha = qcLineAlpha,
                   color = qcCutoffColor,
                   linetype = qcLineType,
                   size = qcLineSize,
                   yintercept = min) +
        scale_y_log10() +
        scale_fill_viridis(discrete = TRUE) +
        theme(axis.text.x = element_text(angle = 90L, hjust = 1L))
    if (isTRUE(metadata(object)[["multiplexedFASTQ"]])) {
        p <- p + facet_wrap(~fileName)
    }
    p
}



.plotUMIsPerCellHistogram <- function(object, min) {
    metrics <- metrics(object)
    p <- ggplot(metrics,
        aes_(x = ~nUMI,
             fill = ~sampleName)) +
        labs(x = "umis per cell") +
        geom_histogram(bins = bins) +
        geom_vline(alpha = qcLineAlpha,
                   color = qcCutoffColor,
                   linetype = qcLineType,
                   size = qcLineSize,
                   xintercept = min) +
        scale_x_log10() +
        scale_y_sqrt() +
        scale_fill_viridis(discrete = TRUE) +
        theme(axis.text.x = element_text(angle = 90L, hjust = 1L))
    if (isTRUE(metadata(object)[["multiplexedFASTQ"]])) {
        p <- p + facet_wrap(~fileName)
    }
    p
}



.plotUMIsPerCell <- function(object, min) {
    plot_grid(.plotUMIsPerCellHistogram(object, min),
              .plotUMIsPerCellBoxplot(object, min),
              labels = "auto",
              nrow = 2L)
}



# Methods ====
#' @rdname plotUMIsPerCell
#' @export
setMethod(
    "plotUMIsPerCell",
    "bcbioSCDataSet",
    function(object, min = 1000L) {
        .plotUMIsPerCell(object, min)
    })



#' @rdname plotUMIsPerCell
#' @export
setMethod(
    "plotUMIsPerCell",
    "bcbioSCFiltered",
    function(object) {
        min <- object %>%
            metadata %>%
            .[["filterParams"]] %>%
            .[["minUMIs"]]
        .plotUMIsPerCell(object, min)
    })
