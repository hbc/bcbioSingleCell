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
#' @param geom Plot type. Supported formats: proportional `histogram`
#'   (*recommended*), raw `ridgeline`, and raw `violin`.
#'
#' @note Here by "cell" we mean "cellular barcode".
#'
#' @return [ggplot].
NULL



# Constructors ====
#' Raw Cellular Barcodes Tibble
#'
#' Used for `geom` parameter: `ridgeline` and `violin` arguments.
#'
#' @author Michael Steinbaugh
#' @keywords internal
#' @noRd
#'
#' @importFrom dplyr group_by left_join mutate
#' @importFrom rlang !! sym
#'
#' @inheritParams plotReadsPerCell
#'
#' @return [tibble] grouped by `sampleName` containing `log10Count` values.
.rawCBTibble <- function(cellularBarcodes, sampleMetadata) {
    sampleMetadata <- sampleMetadata[, c("sampleID", "sampleName")]
    cellularBarcodes %>%
        mutate(log10Count = log10(.data[["nCount"]]),
               cellularBarcode = NULL,
               nCount = NULL) %>%
        left_join(sampleMetadata, by = "sampleID") %>%
        group_by(!!sym("sampleName"))
}



#' Proportional Cellular Barcodes Tibble
#'
#' @author Rory Kirchner, Michael Steinbaugh
#' @keywords internal
#' @noRd
#'
#' @importFrom dplyr bind_rows left_join mutate
#' @importFrom parallel mclapply
#' @importFrom tibble tibble
#'
#' @param rawTibble [.rawCBTibble()] return.
#' @param sampleMetadata [sampleMetadata()] return with `sampleName` columns
#'   that match the `rawTibble`.
#'
#' @details Modified version of Allon Klein Lab MATLAB code.
#'
#' @return [tibble].
.proportionalCBTibble <- function(rawTibble, sampleMetadata) {
    # Ensure `sampleName` is set as factor across both data frames
    rawTibble[["sampleName"]] <-
        as.factor(rawTibble[["sampleName"]])
    sampleMetadata[["sampleName"]] <-
        as.factor(sampleMetadata[["sampleName"]])
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
        mutate(sampleName = as.factor(.data[["sampleName"]])) %>%
        left_join(sampleMetadata, by = "sampleName")
}



#' Plot Raw Cellular Barcodes Violin
#'
#' @author Michael Steinbaugh
#' @keywords internal
#' @noRd
#'
#' @importFrom viridis scale_fill_viridis
#'
#' @inheritParams plotReadsPerCell
#'
#' @return [ggplot].
.plotRawCBViolin <- function(
    tibble,
    interestingGroups = "sampleName",
    cutoffLine = 0,
    multiplexed = FALSE) {
    # Only plot a minimum of 100 reads per cell (2 on X axis). Otherwise the
    # plot gets dominated by cellular barcodes with low read counts.
    tibble <- tibble[tibble[["log10Count"]] >= 2, , drop = FALSE]
    p <- ggplot(
        tibble,
        mapping = aes_string(
            x = "sampleName",
            y = "log10Count",
            fill = interestingGroups)
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
    if (isTRUE(multiplexed) &
        length(unique(tibble[["description"]])) > 1) {
        facets <- c(facets, "description")
    }
    if (isTRUE(.checkAggregate(tibble))) {
        facets <- c(facets, "sampleNameAggregate")
        # Turn off the legend
        p <- p + theme(legend.position = "none")
    }
    if (!is.null(facets)) {
        # Use `free_y` here because of `coord_flip()`
        p <- p + facet_wrap(facets = facets, scales = "free_y")
    }

    p
}



#' Plot Raw Cellular Barcodes Ridgeline
#'
#' @author Michael Steinbaugh
#' @keywords internal
#' @noRd
#'
#' @importFrom ggridges geom_density_ridges
#'
#' @inheritParams plotReadsPerCell
#'
#' @return [ggplot].
.plotRawCBRidgeline <- function(
    tibble,
    interestingGroups = "sampleName",
    cutoffLine = 2,
    multiplexed = FALSE) {
    # Only plot a minimum of 100 reads per cell (2 on X axis). Otherwise the
    # plot gets dominated by cellular barcodes with low read counts.
    tibble <- tibble %>%
        .[.[["log10Count"]] >= 2, , drop = FALSE]
    p <- ggplot(
        tibble,
        mapping = aes_string(
            x = "log10Count",
            y = "sampleName",
            fill = interestingGroups)
    ) +
        labs(title = "raw ridgeline",
             x = "log10 reads per cell") +
        geom_density_ridges(
            alpha = qcPlotAlpha,
            color = lineColor,
            panel_scaling = TRUE,
            scale = 10) +
        scale_fill_viridis(discrete = TRUE) +
        scale_x_sqrt()

    # Cutoff lines
    if (cutoffLine > 0 & length(cutoffLine)) {
        p <- p +
            .qcCutoffLine(xintercept = cutoffLine)
    }

    # Facets
    facets <- NULL
    if (isTRUE(multiplexed) &
        length(unique(tibble[["description"]])) > 1) {
        facets <- c(facets, "description")
    }
    if (isTRUE(.checkAggregate(tibble))) {
        facets <- c(facets, "sampleNameAggregate")
        # Turn off the legend
        p <- p + theme(legend.position = "none")
    }
    if (!is.null(facets)) {
        p <- p + facet_wrap(facets = facets)
    }

    p
}



#' Plot Proportional Cellular Barcodes Histogram
#'
#' @author Michael Steinbaugh
#' @keywords internal
#' @noRd
#'
#' @importFrom viridis scale_color_viridis
#'
#' @inheritParams plotReadsPerCell
#'
#' @return [ggplot].
.plotProportionalCBHistogram <- function(
    tibble,
    interestingGroups = "sampleName",
    cutoffLine = NULL,
    multiplexed = FALSE) {
    p <- ggplot(
        tibble,
        mapping = aes_string(
            x = "log10Count",
            y = "proportion",
            color = interestingGroups)
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
    if (isTRUE(multiplexed) &
        length(unique(tibble[["description"]])) > 1) {
        facets <- c(facets, "description")
    }
    if (isTRUE(.checkAggregate(tibble))) {
        facets <- c(facets, "sampleNameAggregate")
        # Turn off the legend
        p <- p + theme(legend.position = "none")
    }
    if (!is.null(facets)) {
        p <- p + facet_wrap(facets = facets)
    }

    p
}



#' Plot Read Counts per Cell Constructor
#'
#' @author Michael Steinbaugh
#' @keywords internal
#' @noRd
#'
#' @importFrom cowplot draw_plot ggdraw
#'
#' @inherit plotReadsPerCell
.plotReadsPerCell <- function(
    object,
    geom = "histogram",
    interestingGroups,
    filterCells = FALSE) {
    # Currently only supported for `loadSingleCell()` return
    if (metadata(object)[["pipeline"]] != "bcbio") {
        return(warning(paste(
            "'plotReadsPerCell()' currently only supports bcbio pipeline"
        ), call. = FALSE))
    }

    if (isTRUE(filterCells)) {
        object <- .applyFilterCutoffs(object)
    }

    if (missing(interestingGroups)) {
        interestingGroups <- basejump::interestingGroups(object)
    }

    # Acquire the data required for plotting
    cellularBarcodes <- bcbio(object, "cellularBarcodes")
    if (is.null(cellularBarcodes)) {
        stop("Raw cellular barcode counts not saved in object")
    }
    cellularBarcodes <- .bindCellularBarcodes(cellularBarcodes)

    sampleMetadata <- sampleMetadata(
        object,
        interestingGroups = interestingGroups)
    rawTibble <- .rawCBTibble(
        cellularBarcodes = cellularBarcodes,
        sampleMetadata = sampleMetadata)

    if (geom == "histogram") {
        proportionalTibble <- .proportionalCBTibble(
            rawTibble = rawTibble,
            sampleMetadata = sampleMetadata)
    }

    # Need to set the cellular barcode cutoff in log10 to match the plots
    if (!is.null(metadata(object)[["cbCutoff"]])) {
        # This will be deprecated in a future release in favor of the longer
        # spelling variant
        cutoffLine <- metadata(object)[["cbCutoff"]]
    } else if (!is.null(metadata(object)[["cellularBarcodeCutoff"]])) {
        cutoffLine <- metadata(object)[["cellularBarcodeCutoff"]]
    } else {
        warning("Failed to detect cellular barcode cutoff", call. = FALSE)
        cutoffLine <- 0
    }
    cutoffLine <- cutoffLine %>%
        as.numeric() %>%
        log10()

    multiplexed <- metadata(object)[["multiplexedFASTQ"]]

    if (geom == "histogram") {
        p <- .plotProportionalCBHistogram(
            proportionalTibble,
            cutoffLine = cutoffLine,
            multiplexed = multiplexed)
    } else if (geom == "ridgeline") {
        p <- .plotRawCBRidgeline(
            rawTibble,
            cutoffLine = cutoffLine,
            multiplexed = multiplexed)
    } else if (geom == "violin") {
        p <- .plotRawCBViolin(
            rawTibble,
            cutoffLine = cutoffLine,
            multiplexed = multiplexed)
    }

    p
}



# Methods ====
#' @rdname plotReadsPerCell
#' @export
setMethod(
    "plotReadsPerCell",
    signature("bcbioSingleCell"),
    .plotReadsPerCell)
