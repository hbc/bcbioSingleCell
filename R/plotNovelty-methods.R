#' Plot Novelty Score
#'
#' @rdname plotNovelty
#' @name plotNovelty
#' @inherit plotGenesPerCell
#'
#' @details "Novelty" refers to log10 genes detected per count.
NULL



# Constructors ====
.plotNoveltyBoxplot <- function(object, min) {
    metrics <- metrics(object)
    medianNovelty <-
        aggregate(log10GenesPerUMI ~ sampleID, metrics, median) %>%
        left_join(sampleMetadata(object), by = "sampleID")
    interestingGroup <- interestingGroups(object)[[1L]]
    p <- ggplot(metrics,
        aes_(x = ~sampleName,
             y = ~log10GenesPerUMI,
             fill = as.name(interestingGroup))) +
        labs(x = "sample",
             y = "log10 genes per umi") +
        geom_boxplot() +
        geom_hline(alpha = qcLineAlpha,
                   color = qcPassColor,
                   size = qcLineSize,
                   yintercept = min) +
        geom_label(data = medianNovelty,
                   aes_(label = ~round(log10GenesPerUMI, digits = 2L)),
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



.plotNoveltyHistogram <- function(object, min) {
    metrics <- metrics(object)
    p <- ggplot(metrics,
        aes_(x = ~log10GenesPerUMI,
             fill = ~sampleName)) +
        labs(x = "log10 genes per umi") +
        geom_histogram(bins = bins) +
        geom_vline(alpha = qcLineAlpha,
                   color = qcPassColor,
                   size = qcLineSize,
                   xintercept = min) +
        scale_x_sqrt() +
        scale_y_sqrt()
    if (isTRUE(metadata(object)[["multiplexedFASTQ"]])) {
        p <- p + facet_wrap(~fileName)
    }
    p
}



.plotNovelty <- function(object, min) {
    plot_grid(.plotNoveltyHistogram(object, min) +
                  theme(legend.position = "none"),
              .plotNoveltyBoxplot(object, min) +
                  theme(legend.position = "bottom"),
              labels = "auto",
              nrow = 2L)
}



# Methods ====
#' @rdname plotNovelty
#' @export
setMethod(
    "plotNovelty",
    "bcbioSCDataSet",
    function(object, min = 0.8) {
        .plotNovelty(object, min)
    })



#' @rdname plotNovelty
#' @export
setMethod(
    "plotNovelty",
    "bcbioSCSubset",
    function(object) {
        min <- object %>%
            metadata %>%
            .[["filteringCriteria"]] %>%
            .[["minNovelty"]]
        .plotNovelty(object, min)
    })
