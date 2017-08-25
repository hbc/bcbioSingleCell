#' Plot Mitochondrial Transcript Abundance
#'
#' @rdname plotMitoRatio
#' @name plotMitoRatio
NULL



# Constructors ====
.plotMitoRatioBoxplot <- function(object, max) {
    metrics <- metrics(object)
    medianMitoRatio <-
        aggregate(mitoRatio ~ sampleID, metrics, median) %>%
        left_join(sampleMetadata(object), by = "sampleID")
    interestingGroup <- interestingGroups(object)[[1L]]
    p <- ggplot(metrics,
           aes_(x = ~sampleName,
                y = ~mitoRatio * 100L,
                fill = as.name(interestingGroup))) +
        labs(x = "sample",
             y = "% mito counts") +
        geom_boxplot() +
        geom_hline(alpha = qcLineAlpha,
                   color = qcCutoffColor,
                   linetype = qcLineType,
                   size = qcLineSize,
                   yintercept = max * 100L) +
        geom_label(data = medianMitoRatio,
                   aes_(label = ~round(mitoRatio * 100L, digits = 2L)),
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
    if (isTRUE(metadata(object)[["multiplexedFASTQ"]])) {
        p <- p + facet_wrap(~fileName)
    }
    p
}



.plotMitoRatioHistogram <- function(object, max) {
    metrics <- metrics(object)
    p <- ggplot(metrics,
        aes_(x = ~mitoRatio * 100L,
             fill = ~sampleName)) +
        labs(x = "% mito counts") +
        geom_histogram(bins = bins) +
        geom_vline(alpha = qcLineAlpha,
                   color = qcCutoffColor,
                   linetype = qcLineType,
                   size = qcLineSize,
                   xintercept = max * 100L) +
        scale_x_sqrt() +
        scale_y_sqrt() +
        scale_fill_viridis(discrete = TRUE) +
    if (isTRUE(metadata(object)[["multiplexedFASTQ"]])) {
        p <- p + facet_wrap(~fileName)
    }
    p
}



.plotMitoRatioScatterplot <- function(object) {
    metrics <- metrics(object)
    p <- ggplot(metrics,
        aes_(x = ~nCoding,
             y = ~nMito,
             color = ~sampleName)) +
        labs(x = "mito counts",
             y = "coding counts") +
        geom_point(size = 1L) +
        scale_x_sqrt() +
        scale_y_sqrt() +
        scale_color_viridis(discrete = TRUE) +
        theme(axis.text.x = element_text(angle = 90L, hjust = 1L))
    if (isTRUE(metadata(object)[["multiplexedFASTQ"]])) {
        p <- p + facet_wrap(~fileName)
    }
    p
}



.plotMitoRatio <- function(object, max) {
    plot_grid(.plotMitoRatioScatterplot(object),
              .plotMitoRatioHistogram(object, max),
              .plotMitoRatioBoxplot(object, max),
              labels = "auto",
              nrow = 3L)
}



# Methods ====
#' @rdname plotMitoRatio
#' @export
setMethod(
    "plotMitoRatio",
    "bcbioSCDataSet",
    function(object, max = 0.1) {
        .plotMitoRatio(object, max)
    })



#' @rdname plotMitoRatio
#' @export
setMethod(
    "plotMitoRatio",
    "bcbioSCFiltered",
    function(object) {
        max <- object %>%
            metadata %>%
            .[["filterParams"]] %>%
            .[["maxMitoRatio"]]
        .plotMitoRatio(object, max)
    })
