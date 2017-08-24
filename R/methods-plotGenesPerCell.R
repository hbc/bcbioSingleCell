#' Plot Genes per Cell
#'
#' @rdname plotGenesPerCell
#'
#' @param min Recommended minimum value cutoff.
#' @param max Recommended maximum value cutoff.
#'
#' @return [ggplot] grid.



#' @rdname plotGenesPerCell
#' @usage NULL
.plotGenesPerCellBoxplot <- function(object, min, max) {
    metrics <- metrics(object)
    medianGenes <- aggregate(nGene ~ sampleID, metrics, median) %>%
        left_join(sampleMetadata(object), by = "sampleID") %>%
        mutate(nGene = round(.data[["nGene"]]))
    interestingGroup <- interestingGroups(object)[[1L]]
    p <- ggplot(metrics,
                aes_(x = ~sampleName,
                     y = ~nGene,
                     fill = as.name(interestingGroup))) +
        labs(x = "sample",
             y = "genes per cell") +
        geom_boxplot() +
        geom_hline(alpha = qcLineAlpha,
                   color = qcPassColor,
                   size = qcLineSize,
                   yintercept = min) +
        geom_label(data = medianGenes,
                   aes_(label = ~nGene),
                   alpha = qcLabelAlpha,
                   label.padding = unit(0.1, "lines"),
                   show.legend = FALSE) +
        scale_y_sqrt() +
        theme(axis.text.x = element_text(angle = 90L, hjust = 1L))
    if (!is.null(max)) {
        p <- p +
            geom_hline(alpha = qcLineAlpha,
                       color = qcPassColor,
                       size = qcLineSize,
                       yintercept = max)
    }
    if (isTRUE(metadata(object)[["multiplexedFASTQ"]])) {
        p <- p + facet_wrap(~fileName)
    }
    p
}



#' @rdname plotGenesPerCell
#' @usage NULL
.plotGenesPerCellHistogram <- function(object, min, max) {
    metrics <- metrics(object)
    p <- ggplot(metrics,
                aes_(x = ~nGene,
                     fill = ~sampleName)) +
        labs(x = "genes per cell") +
        geom_histogram(bins = bins) +
        geom_vline(alpha = qcLineAlpha,
                   color = qcPassColor,
                   size = qcLineSize,
                   xintercept = min) +
        scale_x_sqrt() +
        scale_y_sqrt() +
        theme(axis.text.x = element_text(angle = 90L, hjust = 1L))
    if (!is.null(max)) {
        p <- p +
            geom_vline(alpha = qcLineAlpha,
                       color = qcPassColor,
                       size = qcLineSize,
                       xintercept = max)
    }
    if (isTRUE(metadata(object)[["multiplexedFASTQ"]])) {
        p <- p + facet_wrap(~fileName)
    }
    p
}



#' @rdname plotGenesPerCell
#' @usage NULL
.plotGenesPerCell <- function(object, min, max) {
    plot_grid(.plotGenesPerCellHistogram(object, min, max) +
                  theme(legend.position = "none"),
              .plotGenesPerCellBoxplot(object, min, max),
              labels = "auto",
              nrow = 2L)
}



#' @rdname plotGenesPerCell
#' @export
setMethod(
    "plotGenesPerCell",
    "bcbioSCDataSet",
    function(object, min = 500L, max = NULL) {
        .plotGenesPerCell(object, min, max)
    })



#' @rdname plotGenesPerCell
#' @export
setMethod(
    "plotGenesPerCell",
    "bcbioSCFiltered",
    function(object) {
        min <- object %>%
            metadata %>%
            .[["filterParams"]] %>%
            .[["minGenes"]]
        max <- object %>%
            metadata %>%
            .[["filterParams"]] %>%
            .[["maxGenes"]]
        .plotGenesPerCell(object, min, max)
    })
