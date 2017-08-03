#' Plot Cellular Barcode Distributions per Sample
#'
#' @rdname plotCellularBarcodes
#' @name plotCellularBarcodes
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @details A violin plot is a comact display of a continuous distribution. It
#'   is a blend of [geom_boxplot()] and [geom_density()]: a violin plot is a
#'   mirrored density plot displayed in the same way as a boxplot.
#'
#' @return [ggplot].
NULL



# Constructors ====
.cbCutoffLine <- function(object) {
    metadata(object)[["cbCutoff"]] %>%
        as.numeric %>%
        log10
}



.plotCBTbl <- function(object) {
    lst <- bcbio(object, "cellularBarcodes")
    if (is.null(cellularBarcodes)) {
        stop("Raw cellular barcode counts not saved in object")
    }
    interestingGroup <- interestingGroups(object)[[1L]]
    meta <- sampleMetadata(object) %>%
        .[, unique(c("fileName",
                             "sampleID",
                             "sampleName",
                             interestingGroup))]
    lst %>%
        .bindCB %>%
        mutate(log10Count = log10(.data[["reads"]]),
               reads = NULL) %>%
        # Only plot barcodes with at least 100 read counts (log10 = 2)
        filter(.data[["log10Count"]] > 2L) %>%
        left_join(meta, by = "sampleID")
}



.plotCBRawViolin <- function(
    plotCBTbl,
    cbCutoffLine,
    multiplexedFASTQ) {
    p <- ggplot(plotCBTbl,
                aes_(x = ~sampleName,
                     y = ~log10Count,
                     fill = ~sampleName)) +
        geom_violin(scale = "width") +
        geom_hline(color = "black",
                   size = qcLineSize,
                   yintercept = cbCutoffLine) +
        labs(title = "raw violin",
             y = "log10 read counts per cell") +
        coord_flip()
    if (isTRUE(multiplexedFASTQ)) {
        p <- p + facet_wrap(~fileName)
    }
    p
}



.plotCBRawHistogram <- function(
    plotCBTbl,
    cbCutoffLine,
    multiplexedFASTQ) {
    p <- ggplot(plotCBTbl,
                aes_(x = ~log10Count,
                     fill = ~sampleName)) +
        labs(title = "raw histogram",
             x = "log10 read counts per cell") +
        geom_histogram(bins = bins) +
        scale_y_sqrt() +
        geom_vline(color = "black",
                   size = qcLineSize,
                   xintercept = cbCutoffLine)
    if (isTRUE(multiplexedFASTQ)) {
        p <- p + facet_wrap(~fileName)
    }
    p
}



.plotCBProportionalHistogram <- function(object) {
    tbl <- .proportionalCB(object)
    cbCutoffLine <- .cbCutoffLine(object)
    p <- ggplot(tbl,
                aes_(x = ~log10Count,
                     y = ~proportion * 100L,
                     color = ~sampleName)) +
        geom_line(alpha = 0.9,
                  size = 1.5) +
        geom_vline(color = "black",
                   size = qcLineSize,
                   xintercept = cbCutoffLine) +
        labs(title = "proportional histogram",
             x = "log10 read counts per cell",
             y = "% of cells")
    if (isTRUE(metadata(object)[["multiplexedFASTQ"]])) {
        p <- p + facet_wrap(~fileName)
    }
    p
}



# Methods ====
#' @rdname plotCellularBarcodes
#' @export
setMethod("plotCellularBarcodes", "bcbioSCDataSet", function(object) {
    tbl <- .plotCBTbl(object)
    cbCutoffLine <- .cbCutoffLine(object)
    multiplexedFASTQ <- metadata(object)[["multiplexedFASTQ"]]

    violin <- .plotCBRawViolin(tbl, cbCutoffLine, multiplexedFASTQ)
    rawHisto <- .plotCBRawHistogram(tbl, cbCutoffLine, multiplexedFASTQ)
    propHisto <- .plotCBProportionalHistogram(object)

    ggdraw() +
        # Coordinates are relative to lower left corner
        draw_plot(violin +
                      xlab("") +
                      theme(legend.position = "none"),
                  x = 0L, y = 0.7, width = 0.5, height = 0.3) +
        draw_plot(rawHisto +
                      ylab("") +
                      theme(axis.text.y = element_blank(),
                            axis.ticks.y = element_blank(),
                            legend.position = "none"),
                  x = 0.5, y = 0.7, width = 0.5, height = 0.3) +
        draw_plot(propHisto +
                      theme(legend.position = "bottom"),
                  x = 0L, y = 0L, width = 1L, height = 0.7)
})
