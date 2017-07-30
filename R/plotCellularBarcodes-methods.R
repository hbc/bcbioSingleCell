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
    cellularBarcodes <- bcbio(object, "cellularBarcodes")
    if (is.null(cellularBarcodes)) {
        stop("Cellular barcode reads not saved in object")
    }
    interestingGroup <- interestingGroups(object)[[1L]]
    meta <- sampleMetadata(object) %>%
        tidy_select(unique(c("fileName",
                             "sampleID",
                             "sampleName",
                             interestingGroup)))
    cellularBarcodes %>%
        .bindCB %>%
        mutate(log10Reads = log10(.data[["reads"]]),
               reads = NULL) %>%
        # Only plot barcodes with at least 100 reads (log10 = 2)
        filter(.data[["log10Reads"]] > 2L) %>%
        left_join(meta, by = "sampleID")
}



.plotCBRawViolin <- function(
    plotCBTbl,
    cbCutoffLine,
    multiplexedFASTQ) {
    p <- ggplot(plotCBTbl,
                aes_(x = ~sampleName,
                     y = ~log10Reads,
                     fill = ~sampleName)) +
        geom_violin(scale = "width") +
        geom_hline(alpha = qcLineAlpha,
                   color = qcFailColor,
                   size = qcLineSize,
                   yintercept = 2L) +
        geom_hline(alpha = qcLineAlpha,
                   color = qcPassColor,
                   size = qcLineSize,
                   yintercept = cbCutoffLine) +
        labs(title = "raw violin",
             y = "log10 reads per cell") +
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
                aes_(x = ~log10Reads,
                     fill = ~sampleName)) +
        labs(title = "raw histogram",
             x = "log10 reads per cell") +
        geom_histogram(bins = bins) +
        scale_y_sqrt() +
        geom_vline(alpha = qcLineAlpha,
                   color = qcFailColor,
                   size = qcLineSize,
                   xintercept = 2L) +
        geom_vline(alpha = qcLineAlpha,
                   color = qcPassColor,
                   size = qcLineSize,
                   xintercept = cbCutoffLine)
    if (isTRUE(multiplexedFASTQ)) {
        p <- p + facet_wrap(~fileName)
    }
    p
}



.plotCBProportionalHistogram <- function(object) {
    cbCutoffLine <- .cbCutoffLine(object)
    p <- ggplot(.proportionalCB(object),
                aes_(x = ~log10ReadsPerCell,
                     y = ~proportionOfCells * 100L,
                     color = ~sampleName)) +
        geom_line() +
        geom_vline(alpha = qcLineAlpha,
                   color = qcFailColor,
                   size = qcLineSize,
                   xintercept = 2L) +
        geom_vline(alpha = qcLineAlpha,
                   color = qcPassColor,
                   size = qcLineSize,
                   xintercept = cbCutoffLine) +
        labs(title = "proportional histogram",
             x = "log10 reads per cell",
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
    # Use defined plotCBTbl here for improved speed
    plotCBTbl <- .plotCBTbl(object)
    cbCutoffLine <- .cbCutoffLine(object)
    multiplexedFASTQ <- metadata(object)[["multiplexedFASTQ"]]
    ggdraw() +
        # Coordinates are relative to lower left corner
        draw_plot(
            .plotCBRawViolin(plotCBTbl,
                             cbCutoffLine,
                             multiplexedFASTQ) +
                xlab("") +
                theme(legend.position = "none"),
            x = 0L, y = 0.7, width = 0.5, height = 0.3) +
        draw_plot(
            .plotCBRawHistogram(plotCBTbl,
                                cbCutoffLine,
                                multiplexedFASTQ) +
                ylab("") +
                theme(axis.text.y = element_blank(),
                      axis.ticks.y = element_blank(),
                      legend.position = "none"),
            x = 0.5, y = 0.7, width = 0.5, height = 0.3) +
        draw_plot(
            .plotCBProportionalHistogram(object) +
                theme(legend.position = "bottom"),
            x = 0L, y = 0L, width = 1L, height = 0.7)
})
