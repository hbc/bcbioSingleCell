#' Plot Read Counts per Cell
#'
#' @name plotReadsPerCell
#' @family Quality Control Functions
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @inheritParams general
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
#' plotReadsPerCell(pbmc_small)
NULL



# Constructors =================================================================
#' Raw Cellular Barcodes Tibble
#'
#' @author Michael Steinbaugh
#' @keywords internal
#' @noRd
#'
#' @return `tibble` grouped by `sampleID`, containing `log10Count` values.
.rawCBTibble <- function(cellularBarcodes, sampleData) {
    sampleData <- sampleData[, c("sampleID", "sampleName")]
    cellularBarcodes %>%
        mutate(
            log10Count = log10(.data[["nCount"]]),
            cellularBarcode = NULL,
            nCount = NULL
        ) %>%
        left_join(sampleData, by = "sampleID") %>%
        group_by(!!sym("sampleID"))
}



#' Proportional Cellular Barcodes Tibble
#'
#' @author Rory Kirchner, Michael Steinbaugh
#' @keywords internal
#' @noRd
#'
#' @param rawTibble [.rawCBTibble()] return.
#' @param sampleData [sampleData()] return with `sampleID` columns
#'   that match the `rawTibble`.
#'
#' @details Modified version of Allon Klein Lab MATLAB code.
#'
#' @return `tibble`.
.proportionalCBTibble <- function(rawTibble, sampleData) {
    # Ensure `sampleID` is set as factor across both data frames
    rawTibble[["sampleID"]] <- as.factor(rawTibble[["sampleID"]])
    sampleData[["sampleID"]] <- as.factor(sampleData[["sampleID"]])
    sampleIDs <- levels(rawTibble[["sampleID"]])
    mclapply(sampleIDs, function(sampleID) {
        cb <- rawTibble %>%
            .[.[["sampleID"]] == sampleID, , drop = FALSE]
        cbHist <- hist(cb[["log10Count"]], n = 100L, plot = FALSE)
        # `counts = fLog` in MATLAB version
        counts <- cbHist[["counts"]]
        # `mids = xLog` in MATLAB version
        mids <-  cbHist[["mids"]]
        tibble(
            "sampleID" = sampleID,
            # log10 reads per cell
            "log10Count" = mids,
            # Proportion of cells
            "proportion" = counts * (10L ^ mids) / sum(counts * (10L ^ mids))
        )
    }) %>%
        bind_rows() %>%
        mutate(sampleID = as.factor(.data[["sampleID"]])) %>%
        left_join(sampleData, by = "sampleID")
}



.plotRawCBViolin <- function(
    data,
    interestingGroups = "sampleName",
    cutoffLine = 2L,
    fill = scale_fill_viridis(discrete = TRUE)
) {
    assertFormalInterestingGroups(data, interestingGroups)
    assertIsANumberOrNULL(cutoffLine)
    assertIsFillScaleDiscreteOrNULL(fill)

    # Only plot a minimum of 100 reads per cell (2 on X axis). Otherwise the
    # plot gets dominated by cellular barcodes with low read counts.
    data <- data[data[["log10Count"]] >= 2L, , drop = FALSE]

    p <- ggplot(
        data = data,
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
            scale = "count"
        ) +
        .medianLabels(data, medianCol = "log10Count", digits = 2L) +
        theme(axis.text.x = element_text(angle = 90L, hjust = 1L))

    # Cutoff lines
    if (cutoffLine > 0L && length(cutoffLine)) {
        p <- p +
            .qcCutoffLine(yintercept = cutoffLine)
    }

    # Color palette
    if (is(fill, "ScaleDiscrete")) {
        p <- p + fill
    }

    # Facets
    facets <- NULL
    if (.isAggregate(data)) {
        facets <- c(facets, "sampleNameAggregate")
    }
    if (is.character(facets)) {
        p <- p + facet_wrap(facets = facets, scales = "free")
    }

    p
}



.plotRawCBRidgeline <- function(
    data,
    interestingGroups = "sampleName",
    cutoffLine = 2L,
    fill = scale_fill_viridis(discrete = TRUE)
) {
    assertFormalInterestingGroups(data, interestingGroups)
    assertIsANumberOrNULL(cutoffLine)
    assertIsFillScaleDiscreteOrNULL(fill)

    # Only plot a minimum of 100 reads per cell (2 on X axis). Otherwise the
    # plot gets dominated by cellular barcodes with low read counts.
    data <- data %>%
        .[.[["log10Count"]] >= 2L, , drop = FALSE]

    p <- ggplot(
        data = data,
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
        .medianLabels(data, medianCol = "log10Count", digits = 2L) +
        theme(axis.text.x = element_text(angle = 90L, hjust = 1L))

    # Cutoff lines
    if (cutoffLine > 0L && length(cutoffLine)) {
        p <- p +
            .qcCutoffLine(xintercept = cutoffLine)
    }

    # Color palette
    if (is(fill, "ScaleDiscrete")) {
        p <- p + fill
    }

    # Facets
    facets <- NULL
    if (.isAggregate(data)) {
        facets <- c(facets, "sampleNameAggregate")
    }
    if (is.character(facets)) {
        p <- p + facet_wrap(facets = facets, scales = "free")
    }

    p
}



.plotProportionalCBHistogram <- function(
    data,
    interestingGroups = "sampleName",
    cutoffLine = NULL,
    color = scale_color_viridis(discrete = TRUE)
) {
    assertFormalInterestingGroups(data, interestingGroups)
    assertIsANumberOrNULL(cutoffLine)
    assertIsColorScaleDiscreteOrNULL(color)

    p <- ggplot(
        data = data,
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
        )

    # Cutoff lines
    if (is.numeric(cutoffLine)) {
        p <- p + .qcCutoffLine(xintercept = cutoffLine)
    }

    # Color palette
    if (is(color, "ScaleDiscrete")) {
        p <- p + color
    }

    # Facets
    facets <- NULL
    if (.isAggregate(data)) {
        facets <- c(facets, "sampleNameAggregate")
    }
    if (is.character(facets)) {
        p <- p + facet_wrap(facets = facets, scales = "free")
    }

    p
}



# Methods ======================================================================
#' @rdname plotReadsPerCell
#' @export
setMethod(
    "plotReadsPerCell",
    signature("bcbioSingleCell"),
    function(
        object,
        geom = c("histogram", "violin", "ridgeline"),
        interestingGroups,
        color = scale_color_viridis(discrete = TRUE),
        fill = scale_fill_viridis(discrete = TRUE)
    ) {
        validObject(object)
        geom <- match.arg(geom)

        if (missing(interestingGroups)) {
            interestingGroups <- bcbioBase::interestingGroups(object)
        }

        # Obtain the sample metadata
        sampleData <- sampleData(object, interestingGroups = interestingGroups)

        # Obtain the cellular barcode distributions
        cbList <- metadata(object)[["cellularBarcodes"]]
        assert_is_list(cbList)

        # Transform cellular barcodes to log10 scale and set up grouping
        rawTibble <- .rawCBTibble(
            cellularBarcodes = .bindCellularBarcodes(cbList),
            sampleData = sampleData
        )
        if (geom == "histogram") {
            proportionalTibble <- .proportionalCBTibble(
                rawTibble = rawTibble,
                sampleData = sampleData
            )
        }

        # Define the log10 cutoff line (to match plot axis)
        cutoffLine <- metadata(object)[["cellularBarcodeCutoff"]]
        assert_is_a_number(cutoffLine)
        cutoffLine <- log10(cutoffLine)

        if (geom == "histogram") {
            .plotProportionalCBHistogram(
                data = proportionalTibble,
                interestingGroups = interestingGroups,
                cutoffLine = cutoffLine,
                color = color
            )
        } else if (geom == "ridgeline") {
            .plotRawCBRidgeline(
                data = rawTibble,
                interestingGroups = interestingGroups,
                cutoffLine = cutoffLine,
                fill = fill
            )
        } else if (geom == "violin") {
            .plotRawCBViolin(
                data = rawTibble,
                interestingGroups = interestingGroups,
                cutoffLine = cutoffLine,
                fill = fill
            )
        }
    }
)



#' @rdname plotReadsPerCell
#' @export
setMethod(
    "plotReadsPerCell",
    signature("seurat"),
    function(object, ...) {
        message("Raw reads not stored in object")
        invisible()
    }
)
