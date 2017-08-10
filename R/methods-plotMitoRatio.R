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
                   color = qcPassColor,
                   size = qcLineSize,
                   yintercept = max * 100L) +
        geom_label(data = medianMitoRatio,
                   aes_(label = ~round(mitoRatio * 100L, digits = 2L)),
                   alpha = qcLabelAlpha,
                   label.padding = unit(0.1, "lines"),
                   show.legend = FALSE) +
        scale_y_sqrt() +
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
                   color = qcPassColor,
                   size = qcLineSize,
                   xintercept = max * 100L) +
        scale_x_sqrt() +
        scale_y_sqrt()
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
        theme(axis.text.x = element_text(angle = 90L, hjust = 1L))
    if (isTRUE(metadata(object)[["multiplexedFASTQ"]])) {
        p <- p + facet_wrap(~fileName)
    }
    p
}



.plotMitoRatio <- function(object, max) {
    plot_grid(.plotMitoRatioScatterplot(object) +
                  theme(legend.position = "none"),
              .plotMitoRatioHistogram(object, max) +
                  theme(legend.position = "none"),
              .plotMitoRatioBoxplot(object, max) +
                  theme(legend.position = "bottom"),
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
    "bcbioSCSubset",
    function(object) {
        max <- object %>%
            metadata %>%
            .[["filterParams"]] %>%
            .[["maxMitoRatio"]]
        .plotMitoRatio(object, max)
    })
