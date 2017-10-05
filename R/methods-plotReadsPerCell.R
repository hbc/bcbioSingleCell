#' Plot Read Counts per Cell
#'
#' @rdname plotReadsPerCell
#' @name plotReadsPerCell
#' @family Quality Control Metrics
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @inheritParams AllGenerics
#' @inheritParams plotGenesPerCell
#'
#' @details A violin plot is a comact display of a continuous distribution. It
#'   is a blend of [geom_boxplot()] and [geom_density()]: a violin plot is a
#'   mirrored density plot displayed in the same way as a boxplot.
#'
#' @note Here by "cell" we mean "cellular barcode".
#'
#' @return [ggplot].
NULL



# Constructors ====
#' Raw Cellular Barcodes
#'
#' @author Michael Steinbaugh
#' @keywords internal
#' @noRd
#'
#' @inheritParams plotReadsPerCell
#'
#' @return [tibble].
.rawCBTibble <- function(
    object,
    filterCells = FALSE) {
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



#' Plot Raw Cellular Barcodes Violin
#'
#' @author Michael Steinbaugh
#' @keywords internal
#' @noRd
#'
#' @inheritParams plotReadsPerCell
#'
#' @return [ggplot].
.plotRawCBViolin <- function(
    tibble,
    cutoffLine = NULL,
    multiplexedFASTQ = FALSE,
    aggregateReplicates = TRUE) {
    # FIXME Is there a way to simplify this in the tibble function instead?
    if (isTRUE(aggregateReplicates) &
        "sampleNameAggregate" %in% colnames(tibble)) {
        xCol <- "sampleNameAggregate"
    } else {
        xCol <- "sampleName"
    }

    p <- ggplot(
        tibble,
        mapping = aes_string(
            x = xCol,
            y = "log10Count",
            fill = "sampleName")
    ) +
        labs(title = "raw violin",
             y = "log10 reads per cell",
             fill = "sample") +
        geom_violin(color = "NA",
                    scale = "width") +
        coord_flip() +
        scale_fill_viridis(discrete = TRUE)

    # Cutoff lines
    if (cutoffLine > 0 & length(cutoffLine)) {
        p <- p +
            .qcCutoffLine(yintercept = cutoffLine)
    }

    # Facets
    facets <- NULL
    if (isTRUE(multiplexedFASTQ)) {
        facets <- c(facets, "fileName")
    }
    if (!isTRUE(aggregateReplicates) &
        "sampleNameAggregate" %in% colnames(tibble)) {
        facets <- c(facets, "sampleNameAggregate")
        # Turn off the legend
        p <- p +
            theme(legend.position = "none")
    }
    if (!is.null(facets)) {
        p <- p +
            # Use `free_y` here because of `coord_flip()`
            facet_wrap(facets = facets, scales = "free_y")
    }

    p
}



#' Plot Raw Cellular Barcodes Histogram
#'
#' @author Michael Steinbaugh
#' @keywords internal
#' @noRd
#'
#' @inheritParams plotReadsPerCell
#'
#' @return [ggplot].
.plotRawCBHistogram <- function(
    tibble,
    cutoffLine = 0,
    multiplexedFASTQ = FALSE,
    aggregateReplicates = TRUE) {
    # FIXME Again, try to set this at the tibble step instead
    if (isTRUE(aggregateReplicates) &
        "sampleNameAggregate" %in% colnames(tibble)) {
        fill <- "sampleNameAggregate"
    } else {
        fill <- "sampleName"
    }

    p <- ggplot(
        tibble,
        mapping = aes_string(
            x = "log10Count",
            fill = fill)
    ) +
        labs(title = "raw histogram",
             x = "log10 reads per cell",
             fill = "sample") +
        geom_histogram(bins = bins) +
        scale_y_sqrt() +
        scale_fill_viridis(discrete = TRUE)

    # Cutoff lines
    if (cutoffLine > 0 & length(cutoffLine)) {
        p <- p +
            .qcCutoffLine(xintercept = cutoffLine)
    }

    # Facets
    facets <- NULL
    if (isTRUE(multiplexedFASTQ)) {
        facets <- c(facets, "fileName")
    }
    if (!isTRUE(aggregateReplicates) &
        "sampleNameAggregate" %in% colnames(tibble)) {
        # FIXME Improve the color contrast here
        facets <- c(facets, "sampleNameAggregate")
        # Turn off the legend
        p <- p +
            theme(legend.position = "none")
    }
    if (!is.null(facets)) {
        p <- p +
            facet_wrap(facets = facets)
    }

    p
}



#' Proportional Cellular Barcodes
#'
#' @author Rory Kirchner, Michael Steinbaugh
#' @keywords internal
#' @noRd
#'
#' @inheritParams plotReadsPerCell
#'
#' @details Modified version of Allon Klein Lab MATLAB code.
#'
#' @return [tibble].
.proportionalCBTibble <- function(
    object,
    filterCells = FALSE,
    aggregateReplicates = TRUE) {
    # FIXME Need to add aggregation support for the calculations here
    if (isTRUE(aggregateReplicates)) {
        stop("aggregateReplicates support not added yet")
    }

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



#' Plot Proportional Cellular Barcodes Histogram
#'
#' @author Michael Steinbaugh
#' @keywords internal
#' @noRd
#'
#' @inheritParams plotReadsPerCell
#'
#' @return [ggplot].
.plotProportionalCBHistogram <- function(
    tibble,
    cutoffLine = NULL,
    multiplexedFASTQ = FALSE,
    aggregateReplicates = TRUE) {

    p <- ggplot(
        tibble,
        mapping = aes_string(
            x = "log10Count",
            y = "proportion",
            color = "sampleName")
    ) +
        geom_line(alpha = 0.8,
                  size = 1.5) +
        labs(title = "proportional histogram",
             x = "log10 reads per cell",
             y = "proportion of cells",
             color = "sample") +
        scale_color_viridis(discrete = TRUE)

    # Cutoff lines
    if (!is.null(cutoffLine) & length(cutoffLine)) {
        p <- p +
            .qcCutoffLine(xintercept = cutoffLine)
    }

    # Facets
    facets <- NULL
    if (isTRUE(multiplexedFASTQ)) {
        facets <- c(facets, "fileName")
    }
    if (!isTRUE(aggregateReplicates) &
        "sampleNameAggregate" %in% colnames(tibble)) {
        facets <- c(facets, "sampleNameAggregate")
        # Turn off the legend
        p <- p +
            theme(legend.position = "none")
    }
    if (!is.null(facets)) {
        p <- p +
            facet_wrap(facets = facets)
    }

    p
}



#' Plot Read Counts per Cell Constructor
#'
#' @author Michael Steinbaugh
#' @keywords internal
#' @noRd
#' @inherit plotReadsPerCell
.plotReadsPerCell <- function(
    object,
    interestingGroup,
    filterCells = FALSE,
    aggregateReplicates = TRUE) {
    if (metadata(object)[["pipeline"]] != "bcbio") {
        warning(paste(
            "'plotReadsPerCell()' currently only supports",
            "bcbio pipeline output for 'bcbioSingleCell' class"),
            call. = FALSE)
        return(NULL)
    }
    if (missing(interestingGroup)) {
        interestingGroup <- interestingGroups(object)[[1]]
    }

    rawTibble <- .rawCBTibble(object, filterCells = filterCells)
    proportionalTibble <-
        .proportionalCBTibble(
            object,
            filterCells = filterCells,
            aggregateReplicates = aggregateReplicates)

    # Need to set the cellular barcode cutoff in log10 to match the plots
    if (!is.null(metadata(object)[["cbCutoff"]])) {
        # This will be deprecated in a future release in favor of the longer
        # spelling variant
        # FIXME Add a warning here if we detect cbCutoff?
        cutoffLine <- metadata(object)[["cbCutoff"]]
    } else if (!is.null(metadata(object)[["cellularBarcodeCutoff"]])) {
        cutoffLine <- metadata(object)[["cellularBarcodeCutoff"]]
    } else {
        warning("Failed to detect cellular barcode cutoff")
        cutoffLine <- 0
    }
    cutoffLine <- cutoffLine %>%
        as.numeric() %>%
        log10()

    multiplexedFASTQ <- metadata(object)[["multiplexedFASTQ"]]

    violin <- .plotRawCBViolin(
            rawTibble,
            cutoffLine = cutoffLine,
            multiplexedFASTQ = multiplexedFASTQ,
            aggregateReplicates = aggregateReplicates)
    rawHistogram <- .plotRawCBHistogram(
            rawTibble,
            cutoffLine = cutoffLine,
            multiplexedFASTQ = multiplexedFASTQ,
            aggregateReplicates = aggregateReplicates)
    proportionalHistogram <- .plotProportionalCBHistogram(
            proportionalTibble,
            cutoffLine = cutoffLine,
            multiplexedFASTQ = multiplexedFASTQ,
            aggregateReplicates = aggregateReplicates)

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



# Methods ====
#' @rdname plotReadsPerCell
#' @export
setMethod("plotReadsPerCell", "bcbioSingleCellANY", .plotReadsPerCell)
