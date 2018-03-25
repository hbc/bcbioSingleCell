# FIXME Allow the user to define the color palette

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
#' load(system.file("extdata/bcb_small.rda", package = "bcbioSingleCell"))
#'
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
    if (.isAggregate(tibble)) {
        facets <- c(facets, "sampleNameAggregate")
    }
    if (is.character(facets)) {
        # Use `free_y` here because of `coord_flip()`
        p <- p + facet_wrap(facets = facets, scales = "free_y")
    }

    p
}



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
    if (.isAggregate(tibble)) {
        facets <- c(facets, "sampleNameAggregate")
    }
    if (is.character(facets)) {
        p <- p + facet_wrap(facets = facets)
    }

    p
}



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
    if (.isAggregate(tibble)) {
        facets <- c(facets, "sampleNameAggregate")
    }
    if (is.character(facets)) {
        p <- p + facet_wrap(facets = facets)
    }

    p
}



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
    sampleData <- sampleData(object, interestingGroups = interestingGroups)

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
        .plotProportionalCBHistogram(
            proportionalTibble,
            interestingGroups = interestingGroups,
            cutoffLine = cutoffLine
        )
    } else if (geom == "ridgeline") {
        .plotRawCBRidgeline(
            rawTibble,
            interestingGroups = interestingGroups,
            cutoffLine = cutoffLine
        )
    } else if (geom == "violin") {
        .plotRawCBViolin(
            rawTibble,
            interestingGroups = interestingGroups,
            cutoffLine = cutoffLine
        )
    }
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
    function(object, ...) {
        message("Raw reads not stored in object")
        invisible()
    }
)
