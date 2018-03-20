# TODO Allow the user to define the color palette

#' Plot Read Counts per Cell
#'
#' @name plotReadsPerCell
#' @family Quality Control Functions
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @inheritParams general
#' @inheritParams plotGenesPerCell
#' @param geom Plot type. Supported formats: proportional `histogram`
#'   (*recommended*), raw `ridgeline`, and raw `violin`.
#'
#' @note Here by "cell" we mean "cellular barcode".
#'
#' @return `ggplot`.
#'
#' @examples
#' # bcbioSingleCell ====
#' plotReadsPerCell(bcb_small)
#'
#' # seurat ====
#' plotReadsPerCell(seurat_small)
NULL



# Constructors =================================================================
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
#' @return `tibble` grouped by `sampleName`, containing `log10Count` values.
.rawCBTibble <- function(cellularBarcodes, sampleData) {
    sampleData <- sampleData[, c("sampleID", "sampleName")]
    cellularBarcodes %>%
        mutate(
            log10Count = log10(.data[["nCount"]]),
            cellularBarcode = NULL,
            nCount = NULL
        ) %>%
        left_join(sampleData, by = "sampleID") %>%
        group_by(!!sym("sampleName"))
}



#' Proportional Cellular Barcodes Tibble
#'
#' @author Rory Kirchner, Michael Steinbaugh
#' @keywords internal
#' @noRd
#'
#' @importFrom dplyr bind_rows left_join mutate
#' @importFrom graphics hist
#' @importFrom parallel mclapply
#' @importFrom tibble tibble
#'
#' @param rawTibble [.rawCBTibble()] return.
#' @param sampleData [sampleData()] return with `sampleName` columns
#'   that match the `rawTibble`.
#'
#' @details Modified version of Allon Klein Lab MATLAB code.
#'
#' @return `tibble`.
.proportionalCBTibble <- function(rawTibble, sampleData) {
    # Ensure `sampleName` is set as factor across both data frames
    rawTibble[["sampleName"]] <-
        as.factor(rawTibble[["sampleName"]])
    sampleData[["sampleName"]] <-
        as.factor(sampleData[["sampleName"]])
    mclapply(
        seq_along(levels(rawTibble[["sampleName"]])), function(a) {
            sampleName <- levels(rawTibble[["sampleName"]])[[a]]
            cb <- rawTibble %>%
                .[.[["sampleName"]] == sampleName, , drop = FALSE]
            cbHist <- hist(cb[["log10Count"]], n = 100L, plot = FALSE)
            # `counts = fLog` in MATLAB version
            counts <- cbHist[["counts"]]
            # `mids = xLog` in MATLAB version
            mids <-  cbHist[["mids"]]
            tibble(
                sampleName = sampleName,
                # log10 reads per cell
                log10Count = mids,
                # Proportion of cells
                proportion = counts * (10L ^ mids) / sum(counts * (10L ^ mids))
            )
        }) %>%
        bind_rows() %>%
        mutate(sampleName = as.factor(.data[["sampleName"]])) %>%
        left_join(sampleData, by = "sampleName")
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
#' @return `ggplot`.
.plotRawCBViolin <- function(
    tibble,
    interestingGroups = "sampleName",
    cutoffLine = 0L
) {
    # Only plot a minimum of 100 reads per cell (2 on X axis). Otherwise the
    # plot gets dominated by cellular barcodes with low read counts.
    tibble <- tibble[tibble[["log10Count"]] >= 2L, , drop = FALSE]
    p <- ggplot(
        tibble,
        mapping = aes_string(
            x = "sampleName",
            y = "log10Count",
            fill = interestingGroups
        )
    ) +
        labs(y = "log10 reads per cell") +
        geom_violin(
            alpha = qcPlotAlpha,
            color = lineColor,
            scale = "width"
        ) +
        coord_flip() +
        scale_fill_viridis(discrete = TRUE)

    # Cutoff lines
    if (cutoffLine > 0L && length(cutoffLine)) {
        p <- p +
            .qcCutoffLine(yintercept = cutoffLine)
    }

    # Facets
    facets <- NULL
    if (isTRUE(.checkAggregate(tibble))) {
        facets <- c(facets, "sampleNameAggregate")
    }
    if (is.character(facets)) {
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
#' @return `ggplot`.
.plotRawCBRidgeline <- function(
    tibble,
    interestingGroups = "sampleName",
    cutoffLine = 2L) {
    # Only plot a minimum of 100 reads per cell (2 on X axis). Otherwise the
    # plot gets dominated by cellular barcodes with low read counts.
    tibble <- tibble %>%
        .[.[["log10Count"]] >= 2L, , drop = FALSE]
    p <- ggplot(
        tibble,
        mapping = aes_string(
            x = "log10Count",
            y = "sampleName",
            fill = interestingGroups
        )
    ) +
        labs(x = "log10 reads per cell") +
        geom_density_ridges(
            alpha = qcPlotAlpha,
            color = lineColor,
            panel_scaling = TRUE,
            scale = 10L
        ) +
        scale_fill_viridis(discrete = TRUE) +
        scale_x_sqrt()

    # Cutoff lines
    if (cutoffLine > 0L && length(cutoffLine)) {
        p <- p +
            .qcCutoffLine(xintercept = cutoffLine)
    }

    # Facets
    facets <- NULL
    if (isTRUE(.checkAggregate(tibble))) {
        facets <- c(facets, "sampleNameAggregate")
    }
    if (is.character(facets)) {
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
#' @return `ggplot`.
.plotProportionalCBHistogram <- function(
    tibble,
    interestingGroups = "sampleName",
    cutoffLine = NULL) {
    p <- ggplot(
        tibble,
        mapping = aes_string(
            x = "log10Count",
            y = "proportion",
            color = interestingGroups
        )
    ) +
        geom_line(
            alpha = qcPlotAlpha,
            size = 1.5
        ) +
        labs(
            x = "log10 reads per cell",
            y = "proportion of cells"
        ) +
        scale_color_viridis(discrete = TRUE)

    # Cutoff lines
    if (!is.null(cutoffLine) && length(cutoffLine)) {
        p <- p +
            .qcCutoffLine(xintercept = cutoffLine)
    }

    # Facets
    facets <- NULL
    if (isTRUE(.checkAggregate(tibble))) {
        facets <- c(facets, "sampleNameAggregate")
    }
    if (is.character(facets)) {
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
#' @importFrom bcbioBase interestingGroups
#' @importFrom cowplot draw_plot ggdraw
#'
#' @inherit plotReadsPerCell
.plotReadsPerCell <- function(
    object,
    geom = "histogram",
    interestingGroups
) {
    validObject(object)
    skipMessage <- "Object doesn't contain unfiltered cellular barcodes"

    if (missing(interestingGroups)) {
        interestingGroups <- bcbioBase::interestingGroups(object)
    }

    # Obtain the cellular barcode distributions
    cellularBarcodes <- metadata(object)[["cellularBarcodes"]]
    if (is.null(cellularBarcodes)) {
        return(inform(skipMessage))
    }
    if (is.list(cellularBarcodes)) {
        cellularBarcodes <- .bindCellularBarcodes(cellularBarcodes)
    }

    # Obtain the sample metadata
    sampleData <- sampleData(
        object,
        interestingGroups = interestingGroups
    )

    # Mutate cellular barcodes to log10 and set up grouping
    rawTibble <- .rawCBTibble(
        cellularBarcodes = cellularBarcodes,
        sampleData = sampleData
    )

    if (geom == "histogram") {
        proportionalTibble <- .proportionalCBTibble(
            rawTibble = rawTibble,
            sampleData = sampleData
        )
    }

    # Define the log10 cutoff line (to match plot axis)
    if (!is.null(metadata(object)[["cellularBarcodeCutoff"]])) {
        cutoffLine <- metadata(object)[["cellularBarcodeCutoff"]]
    } else {
        warn("Cellular barcode cutoff is not saved in object")
        cutoffLine <- 0L
    }
    cutoffLine <- cutoffLine %>%
        as.numeric() %>%
        log10()

    if (geom == "histogram") {
        p <- .plotProportionalCBHistogram(
            proportionalTibble,
            interestingGroups = interestingGroups,
            cutoffLine = cutoffLine
        )
    } else if (geom == "ridgeline") {
        p <- .plotRawCBRidgeline(
            rawTibble,
            interestingGroups = interestingGroups,
            cutoffLine = cutoffLine
        )
    } else if (geom == "violin") {
        p <- .plotRawCBViolin(
            rawTibble,
            interestingGroups = interestingGroups,
            cutoffLine = cutoffLine
        )
    }

    p
}



# Methods ======================================================================
#' @rdname plotReadsPerCell
#' @export
setMethod(
    "plotReadsPerCell",
    signature("bcbioSingleCell"),
    .plotReadsPerCell
)



#' @rdname plotReadsPerCell
#' @export
setMethod(
    "plotReadsPerCell",
    signature("seurat"),
    .plotReadsPerCell
)
