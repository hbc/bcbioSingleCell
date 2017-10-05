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
#' @return [tibble] grouped by `sampleName` containing `log10Count` values.
.rawCBTibble <- function(
    object,
    filterCells = FALSE,
    aggregateReplicates = TRUE) {
    cellularBarcodes <- bcbio(object, "cellularBarcodes")
    if (is.null(cellularBarcodes)) {
        stop("Raw cellular barcode counts not saved in object")
    }
    cellularBarcodes <- .bindCellularBarcodes(cellularBarcodes)
    # Keep only the cells that have passed filtering, if desired
    if (isTRUE(filterCells)) {
        cells <- metadata(object)[["filterCells"]]
        if (!is.null(cells)) {
            cellularBarcodes <- cellularBarcodes %>%
                .[.[["cellID"]] %in% cells, , drop = FALSE]
        }
    }
    meta <- metadata(object)[["sampleMetadata"]]
    if (isTRUE(aggregateReplicates) &
        "sampleNameAggregate" %in% colnames(meta)) {
        meta[["sampleName"]] <- meta[["sampleNameAggregate"]]
        meta[["sampleNameAggregate"]] <- NULL
    }
    meta <- meta[, metaPriorityCols]
    cellularBarcodes %>%
        mutate(log10Count = log10(.data[["nCount"]]),
               cellularBarcode = NULL,
               nCount = NULL) %>%
        left_join(meta, by = "sampleID") %>%
        .[, c("sampleName", "log10Count")] %>%
        group_by(!!sym("sampleName"))
}



#' Proportional Cellular Barcodes
#'
#' @author Rory Kirchner, Michael Steinbaugh
#' @keywords internal
#' @noRd
#'
#' @param rawTibble [.rawCBTibble()] return.
#' @param perSampleMetadata [sampleMetadata()] return with `sampleName` columns
#'   that match the `rawTibble`.
#'
#' @details Modified version of Allon Klein Lab MATLAB code.
#'
#' @return [tibble].
.proportionalCBTibble <- function(
    rawTibble,
    perSampleMetadata) {
    mclapply(
        seq_along(levels(rawTibble[["sampleName"]])), function(a) {
            sampleName <- levels(rawTibble[["sampleName"]])[[a]]
            cb <- rawTibble %>%
                .[.[["sampleName"]] == sampleName, , drop = FALSE]
            cbHist <- hist(cb[["log10Count"]], n = 100, plot = FALSE)
            # `counts = fLog` in MATLAB version
            counts <- cbHist[["counts"]]
            # `mids = xLog` in MATLAB version
            mids <-  cbHist[["mids"]]
            tibble(
                sampleName = sampleName,
                # log10 reads per cell
                log10Count = mids,
                # Proportion of cells
                proportion = counts * (10 ^ mids) / sum(counts * (10 ^ mids)))
        }) %>%
        bind_rows() %>%
        mutate(sampleName = factor(.data[["sampleName"]])) %>%
        left_join(perSampleMetadata, by = "sampleName")
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
    interestingGroup = "sampleName",
    cutoffLine = 0,
    multiplexedFASTQ = FALSE,
    aggregateReplicates = TRUE) {
    # Only plot a minimum of 100 reads per cell (2 on X axis). Otherwise the
    # plot gets dominated by cellular barcodes with low read counts.
    tibble <- tibble %>%
        .[.[["log10Count"]] >= 2, , drop = FALSE]
    p <- ggplot(
        tibble,
        mapping = aes_string(
            x = "sampleName",
            y = "log10Count",
            fill = interestingGroup)
    ) +
        labs(title = "raw violin",
             y = "log10 reads per cell") +
        geom_violin(
            alpha = qcPlotAlpha,
            color = lineColor,
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
.plotRawCBRidgeline <- function(
    tibble,
    interestingGroup = "sampleName",
    cutoffLine = 2,
    multiplexedFASTQ = FALSE,
    aggregateReplicates = TRUE) {
    # Only plot a minimum of 100 reads per cell (2 on X axis). Otherwise the
    # plot gets dominated by cellular barcodes with low read counts.
    tibble <- tibble %>%
        .[.[["log10Count"]] >= 2, , drop = FALSE]
    p <- ggplot(
        tibble,
        mapping = aes_string(
            x = "log10Count",
            y = "sampleName",
            fill = interestingGroup)
    ) +
        labs(title = "raw histogram",
             x = "log10 reads per cell") +
        geom_density_ridges(
            alpha = qcPlotAlpha,
            color = lineColor,
            panel_scaling = TRUE,
            scale = qcRidgeScale) +
        scale_fill_viridis(discrete = TRUE) +
        scale_x_sqrt()

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
    interestingGroup = "sampleName",
    cutoffLine = NULL,
    multiplexedFASTQ = FALSE,
    aggregateReplicates = TRUE) {
    p <- ggplot(
        tibble,
        mapping = aes_string(
            x = "log10Count",
            y = "proportion",
            color = interestingGroup)
    ) +
        geom_line(alpha = qcPlotAlpha,
                  size = 1.5) +
        labs(title = "proportional histogram",
             x = "log10 reads per cell",
             y = "proportion of cells") +
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

    rawTibble <- .rawCBTibble(
        object,
        filterCells = filterCells)
    perSampleMetadata <- sampleMetadata(
        object, aggregateReplicates = aggregateReplicates)
    proportionalTibble <- .proportionalCBTibble(
        rawTibble = rawTibble,
        perSampleMetadata = perSampleMetadata)

    # Need to set the cellular barcode cutoff in log10 to match the plots
    if (!is.null(metadata(object)[["cbCutoff"]])) {
        # This will be deprecated in a future release in favor of the longer
        # spelling variant
        warning(paste(
            "'cbCutoff' in metadata will be deprecated in favor of",
            "'cellularBarcodeCutoff' in a future update"
        ))
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
    ridgeline <- .plotRawCBRidgeline(
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
        suppresMessages(draw_plot(
            ridgeline +
                ylab("") +
                theme(axis.text.y = element_blank(),
                      axis.ticks.y = element_blank(),
                      legend.position = "none"),
            x = 0.5, y = 0.7, width = 0.5, height = 0.3)) +
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
