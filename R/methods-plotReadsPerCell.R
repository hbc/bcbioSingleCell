#' Plot Read Counts per Cell
#'
#' @rdname plotReadsPerCell
#' @name plotReadsPerCell
#' @family Quality Control Metrics
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @inheritParams AllGenerics
#' @param interestingGroup Interesting group, used for colors.
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
.rawCBTibble <- function(object, filterCells = FALSE) {
    cellularBarcodes <- bcbio(object, "cellularBarcodes")
    if (is.null(cellularBarcodes)) {
        stop("Raw cellular barcode counts not saved in object")
    }
    cellularBarcodes <- .bindCellularBarcodes(cellularBarcodes)
    # Keep only the cells that have passed filtering, if desired
    if (isTRUE(filterCells)) {
        cells <- metadata(object)[["filterCells"]]
        cellularBarcodes <- cellularBarcodes %>%
            .[.[["cellID"]] %in% cells, , drop = FALSE]
    }
    cellularBarcodes %>%
        mutate(log10Count = log10(.data[["nCount"]]),
               cellularBarcode = NULL,
               nCount = NULL) %>%
        # Only plot barcodes with at least 100 read counts (log10 = 2)
        dplyr::filter(.data[["log10Count"]] > 2) %>%
        left_join(sampleMetadata(object), by = "sampleID")
}



.plotRawCBViolin <- function(
    tibble,
    cutoffLine,
    multiplexedFASTQ) {
    p <- ggplot(
        tibble,
        mapping = aes_string(
            x = "sampleName",
            y = "log10Count",
            fill = "sampleName")
    ) +
        geom_violin(color = "NA", scale = "width") +
        labs(title = "raw violin",
             y = "log10 reads per cell") +
        coord_flip() +
        scale_fill_viridis(discrete = TRUE)
    if (!is.null(cutoffLine) & length(cutoffLine)) {
        p <- p +
            .qcCutoffLine(yintercept = cutoffLine)
    }
    if (isTRUE(multiplexedFASTQ)) {
        p <- p +
            facet_wrap(~fileName)
    }
    p
}



.plotRawCBHistogram <- function(
    tibble,
    cutoffLine,
    multiplexedFASTQ) {
    p <- ggplot(
        tibble,
        mapping = aes_string(
            x = "log10Count",
            fill = "sampleName")
    ) +
        labs(title = "raw histogram",
             x = "log10 reads per cell") +
        geom_histogram(bins = bins) +
        scale_y_sqrt() +
        scale_fill_viridis(discrete = TRUE)
    if (!is.null(cutoffLine) & length(cutoffLine)) {
        p <- p +
            .qcCutoffLine(xintercept = cutoffLine)
    }
    if (isTRUE(multiplexedFASTQ)) {
        p <- p +
            facet_wrap(~fileName)
    }
    p
}



#' Proportional Cellular Barcodes
#'
#' @author Rory Kirchner, Michael Steinbaugh
#'
#' @inheritParams AllGenerics
#' @param filterCells Use only the cells that have passed filtering.
#'
#' @details Modified version of Klein Lab MATLAB code.
#'
#' @return [tibble].
#' @noRd
.proportionalCBTibble <- function(object, filterCells = FALSE) {
    cellularBarcodes <- bcbio(object, "cellularBarcodes")
    lapply(seq_along(cellularBarcodes), function(a) {
        sampleID <- names(cellularBarcodes)[[a]]
        cb <- cellularBarcodes[[a]] %>%
            mutate(
                cellID = paste(sampleID,
                               .data[["cellularBarcode"]],
                               sep = "_"),
                log10Count = log10(.data[["nCount"]]),
                cellularBarcode = NULL,
                nCount = NULL)
        # Keep only the cells that have passed filtering, if desired
        if (isTRUE(filterCells)) {
            cells <- metadata(object)[["filterCells"]]
            cb <- cb %>%
                .[.[["cellID"]] %in% cells, , drop = FALSE]
            # Return NULL if none of the cells in the sample passed filtering
            if (!nrow(cb)) {
                warning(paste("No cells passed filtering for", sampleID),
                        call. = FALSE)
                return(NULL)
            }
        }
        cbHist <- hist(cb[["log10Count"]], n = 100, plot = FALSE)
        # `fLog` in MATLAB version
        counts <- cbHist[["counts"]]
        # `xLog` in MATLAB version
        mids <-  cbHist[["mids"]]
        tibble(
            sampleID = sampleID,
            # log10 reads per cell
            log10Count = mids,
            # Proportion of cells
            proportion = counts * (10 ^ mids) /
                sum(counts * (10 ^ mids)))
    }) %>%
        bind_rows() %>%
        left_join(sampleMetadata(object), by = "sampleID")
}



.plotProportionalCBHistogram <- function(
    tibble,
    cutoffLine,
    multiplexedFASTQ) {
    p <- ggplot(
        tibble,
        mapping = aes_(
            x = ~log10Count,
            y = ~proportion * 100,
            color = ~sampleName)) +
        geom_line(alpha = 0.9,
                  size = 1.5) +
        labs(title = "proportional histogram",
             x = "log10 reads per cell",
             y = "% of cells") +
        scale_color_viridis(discrete = TRUE)
    if (!is.null(cutoffLine) & length(cutoffLine)) {
        p <- p +
            .qcCutoffLine(xintercept = cutoffLine)
    }
    if (isTRUE(multiplexedFASTQ)) {
        p <- p +
            facet_wrap(~fileName)
    }
    p
}



.plotCellularBarcode <- function(
    violin,
    rawHistogram,
    proportionalHistogram) {
    ggdraw() +
        # Coordinates are relative to lower left corner
        draw_plot(
            violin +
                xlab("") +
                theme(legend.position = "none"),
            x = 0, y = 0.7, width = 0.5, height = 0.3) +
        draw_plot(
            rawHistogram +
                ylab("") +
                theme(axis.text.y = element_blank(),
                      axis.ticks.y = element_blank(),
                      legend.position = "none"),
            x = 0.5, y = 0.7, width = 0.5, height = 0.3) +
        draw_plot(
            proportionalHistogram +
                theme(legend.justification = "center",
                      legend.position = "bottom"),
            x = 0, y = 0, width = 1, height = 0.7)
}



.plotReadsPerCell <- function(
    object,
    interestingGroup,
    filterCells = FALSE) {
    # This function currently only supports bcbio pipeline
    if (metadata(object)[["pipeline"]] != "bcbio") {
        warning(paste(
            "'plotReadsPerCell()' currently only supports",
            "bcbio pipeline for 'bcbioSingleCell' class"),
            call. = FALSE)
        return(NULL)
    }

    if (missing(interestingGroup)) {
        interestingGroup <- interestingGroups(object)[[1]]
    }

    rawTibble <-
        .rawCBTibble(object, filterCells = filterCells)
    proportionalTibble <-
        .proportionalCBTibble(object, filterCells = filterCells)

    # Need to set the cellular barcode cutoff in log10 to match the plots
    if (!is.null(metadata(object)[["cbCutoff"]])) {
        # This will be deprecated in a future release in favor of the longer
        # spelling variant
        cutoffLine <- metadata(object)[["cbCutoff"]]
    } else if (!is.null(metadata(object)[["cellularBarcodeCutoff"]])) {
        cutoffLine <- metadata(object)[["cellularBarcodeCutoff"]]
    } else {
        warning("Failed to detect cellular barcode cutoff")
        cutoffLine <- NULL
    }
    cutoffLine <- cutoffLine %>%
        as.numeric() %>%
        log10()

    multiplexedFASTQ <- metadata(object)[["multiplexedFASTQ"]]

    .plotCellularBarcode(
        .plotRawCBViolin(
            rawTibble,
            cutoffLine = cutoffLine,
            multiplexedFASTQ = multiplexedFASTQ),
        .plotRawCBHistogram(
            rawTibble,
            cutoffLine = cutoffLine,
            multiplexedFASTQ = multiplexedFASTQ),
        .plotProportionalCBHistogram(
            proportionalTibble,
            cutoffLine = cutoffLine,
            multiplexedFASTQ = multiplexedFASTQ))
}



# Methods ====
#' @rdname plotReadsPerCell
#' @export
setMethod("plotReadsPerCell", "bcbioSingleCellANY", .plotReadsPerCell)
