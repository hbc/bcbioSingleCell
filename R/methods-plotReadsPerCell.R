#' Plot Read Counts per Cell
#'
#' @rdname plotReadsPerCell
#' @name plotReadsPerCell
#'
#' @details A violin plot is a comact display of a continuous distribution. It
#'   is a blend of [geom_boxplot()] and [geom_density()]: a violin plot is a
#'   mirrored density plot displayed in the same way as a boxplot.
#'
#' @note Here by cell we mean "cellular barcode".
#'
#' @return [ggplot].
NULL



# Constructors ====
.cbCutoffLine <- function(object) {
    metadata(object)[["cbCutoff"]] %>%
        as.numeric %>%
        log10
}



.cbTblFromList <- function(object) {
    lst <- bcbio(object, "cellularBarcodes")
    if (is.null(lst)) {
        stop("Raw cellular barcode counts not saved in object")
    }
    interestingGroup <- interestingGroups(object)[[1L]]
    meta <- sampleMetadata(object) %>%
        .[, unique(c(metaPriorityCols, interestingGroup))]
    lst %>%
        .bindCB %>%
        mutate(log10Count = log10(.data[["nCount"]]),
               nCount = NULL) %>%
        # Only plot barcodes with at least 100 read counts (log10 = 2)
        tidy_filter(.data[["log10Count"]] > 2L) %>%
        left_join(meta, by = "sampleID")
}



.cbTblFromMetrics <- function(object) {
    interestingGroup <- interestingGroups(object)[[1L]]
    meta <- sampleMetadata(object) %>%
        .[, unique(c(metaPriorityCols, interestingGroup))]
    metrics(object) %>%
        as("tibble") %>%
        mutate(log10Count = log10(.data[["nCount"]]),
               nCount = NULL) %>%
        tidy_filter(.data[["log10Count"]] > 2L)
}



.plotCBRawViolin <- function(
    object,
    cutoffLine = NULL,
    multiplexedFASTQ = FALSE) {
    p <- ggplot(object,
                aes_(x = ~sampleName,
                     y = ~log10Count,
                     fill = ~sampleName)) +
        geom_violin(scale = "width") +
        labs(title = "raw violin",
             y = "log10 reads per cell") +
        coord_flip()
    if (!is.null(cutoffLine) & length(cutoffLine)) {
        p <- p + geom_hline(color = "black",
                            size = qcLineSize,
                            yintercept = cutoffLine)
    }
    if (isTRUE(multiplexedFASTQ)) {
        p <- p + facet_wrap(~fileName)
    }
    p
}



.plotCBRawHisto <- function(
    object,
    cutoffLine = NULL,
    multiplexedFASTQ = FALSE) {
    p <- ggplot(object,
                aes_(x = ~log10Count,
                     fill = ~sampleName)) +
        labs(title = "raw histogram",
             x = "log10 reads per cell") +
        geom_histogram(bins = bins) +
        scale_y_sqrt()
    if (!is.null(cutoffLine) & length(cutoffLine)) {
        p <- p + geom_vline(color = "black",
                            size = qcLineSize,
                            xintercept = cutoffLine)
    }
    if (isTRUE(multiplexedFASTQ)) {
        p <- p + facet_wrap(~fileName)
    }
    p
}



#' Proportional Cellular Barcodes
#'
#' @rdname internal-propCB
#' @keywords internal
#'
#' @inheritParams AllGenerics
#'
#' @details Modified version of Klein Lab MATLAB code.
#'
#' @return [tibble].
.propTblFromDataSet <- function(object) {
    metadata <- sampleMetadata(object) %>%
        .[, metaPriorityCols]
    lst <- bcbio(object, "cellularBarcodes")
    lapply(seq_along(lst), function(a) {
        cb <- lst[[a]] %>%
            mutate(log10Count = log10(.data[["nCount"]]))
        cbHist <- hist(cb[["log10Count"]], n = 100L, plot = FALSE)
        # `fLog` in Klein Lab code
        counts <- cbHist[["counts"]]
        # `xLog` in Klein Lab code
        mids <-  cbHist[["mids"]]
        tibble(
            sampleID = names(lst)[[a]],
            # log10 reads per cell
            log10Count = mids,
            # Proportion of cells
            proportion = counts * (10L ^ mids) /
                sum(counts * (10L ^ mids)))
    }) %>%
        set_names(names(lst)) %>%
        bind_rows %>%
        left_join(metadata, by = "sampleID")
}



.propTblFromSubset <- function(object) {
    metadata <- sampleMetadata(object) %>%
        .[, metaPriorityCols]
    metrics <- metrics(object) %>%
        mutate(log10Count = log10(.data[["nCount"]])) %>%
        .[, c(metaPriorityCols, "log10Count")]
    uniques <- unique(metrics[["sampleID"]])
    lapply(seq_along(uniques), function(a) {
        cb <- metrics %>%
            .[.[["sampleID"]] %in% uniques[[a]], , drop = FALSE]
        cbHist <- hist(cb[["log10Count"]], n = 100L, plot = FALSE)
        # `fLog` in Klein Lab code
        counts <- cbHist[["counts"]]
        # `xLog` in Klein Lab code
        mids <-  cbHist[["mids"]]
        tibble(
            sampleID = uniques[[a]],
            # log10 reads per cell
            log10Count = mids,
            # Proportion of cells
            proportion = counts * (10L ^ mids) /
                sum(counts * (10L ^ mids)))
    }) %>%
        set_names(uniques) %>%
        bind_rows %>%
        left_join(metadata, by = "sampleID")
}



.plotCBPropHisto <- function(
    tbl, cutoffLine = NULL, multiplexedFASTQ = FALSE) {
    p <- ggplot(tbl,
                aes_(x = ~log10Count,
                     y = ~proportion * 100L,
                     color = ~sampleName)) +
        geom_line(alpha = 0.9,
                  size = 1.5) +
        labs(title = "proportional histogram",
             x = "log10 reads per cell",
             y = "% of cells")
    if (!is.null(cutoffLine) & length(cutoffLine)) {
        p <- p + geom_vline(color = "black",
                            size = qcLineSize,
                            xintercept = cutoffLine)
    }
    if (isTRUE(multiplexedFASTQ)) {
        p <- p + facet_wrap(~fileName)
    }
    p
}



.plotCB <- function(violin, rawHisto, propHisto) {
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
}



# Methods ====
#' @rdname plotReadsPerCell
#' @export
setMethod("plotReadsPerCell", "bcbioSCDataSet", function(object) {
    # Currently only supports bcbio pipeline
    if (metadata(object)[["pipeline"]] != "bcbio") {
        warning(paste("`plotReadsPerCell()` currently only supports",
                      "bcbio pipeline for `bcbioSCDataSet` class"),
                call. = FALSE)
        return(NULL)
    }
    rawTbl <- .cbTblFromList(object)
    propTbl <- .propTblFromDataSet(object)
    cutoffLine <- .cbCutoffLine(object)
    multiplexedFASTQ <- metadata(object)[["multiplexedFASTQ"]]
    .plotCB(
        .plotCBRawViolin(
            rawTbl,
            cutoffLine = cutoffLine,
            multiplexedFASTQ = multiplexedFASTQ),
        .plotCBRawHisto(
            rawTbl,
            cutoffLine = cutoffLine,
            multiplexedFASTQ = multiplexedFASTQ),
        .plotCBPropHisto(
            propTbl,
            cutoffLine = cutoffLine,
            multiplexedFASTQ = multiplexedFASTQ))
})



#' @rdname plotReadsPerCell
#' @export
setMethod("plotReadsPerCell", "bcbioSCSubset", function(object) {
    rawTbl <- .cbTblFromMetrics(object)
    propTbl <- .propTblFromSubset(object)
    multiplexedFASTQ <- metadata(object)[["multiplexedFASTQ"]]
    .plotCB(
        .plotCBRawViolin(
            rawTbl,
            multiplexedFASTQ = multiplexedFASTQ),
        .plotCBRawHisto(
            rawTbl,
            multiplexedFASTQ = multiplexedFASTQ),
        .plotCBPropHisto(
            propTbl,
            multiplexedFASTQ = multiplexedFASTQ))
})
