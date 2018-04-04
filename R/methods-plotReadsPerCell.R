#' Plot Read Counts per Cell
#'
#' Plot the distribution of read counts for all unfiltered cellular barcodes.
#'
#' @name plotReadsPerCell
#' @family Quality Control Functions
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @inheritParams general
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
.readsPerCell <- function(object, interestingGroups) {
    validObject(object)
    if (missing(interestingGroups)) {
        interestingGroups <- bcbioBase::interestingGroups(object)
    }

    # Obtain the sample metadata
    sampleData <- sampleData(
        object = object,
        interestingGroups = interestingGroups,
        return = "data.frame"
    )

    # Obtain the read counts. Use the unfiltered reads stashed in the metadata
    # if available, otherwise use the metrics return.
    cbList <- metadata(object)[["cellularBarcodes"]]
    if (is.list(cbList)) {
        data <- .bindCellularBarcodes(cbList)
    } else {
        data <- metrics(object)
    }

    # Return `NULL` if `nCount` isn't present
    if (!"nCount" %in% colnames(data)) {
        warn("object does not contain nCount column in `metrics()`")
        return(NULL)
    }

    left_join(
        x = data[, c("sampleID", "nCount")],
        y = sampleData,
        by = "sampleID"
    ) %>%
        as_tibble() %>%
        group_by(!!sym("sampleID"))
}



.plotReadsPerCell <- function(
    object,
    interestingGroups,
    geom = c("histogram", "ecdf", "ridgeline", "violin"),
    color = scale_color_viridis(discrete = TRUE),
    fill = scale_fill_viridis(discrete = TRUE)
) {
    # Passthrough: color, fill
    validObject(object)
    if (missing(interestingGroups)) {
        interestingGroups <- bcbioBase::interestingGroups(object)
    }
    geom <- match.arg(geom)

    data <- .readsPerCell(
        object = object,
        interestingGroups = interestingGroups
    )
    if (!is.data.frame(data)) {
        return(invisible())
    }

    cutoff <- metadata(object)[["cellularBarcodeCutoff"]]
    assert_is_an_integer(cutoff)

    if (geom == "histogram") {
        sampleData <- sampleData(
            object = object,
            interestingGroups = interestingGroups,
            return = "data.frame"
        )
        data <- .proportionalReadsPerCell(
            data = data,
            sampleData = sampleData
        )
        p <- .plotReadsPerCellHistogram(
            data = data,
            cutoff = cutoff,
            color = color
        )
    } else if (geom == "ecdf") {
        p <- .plotReadsPerCellECDF(
            data = data,
            cutoff = cutoff,
            color = color
        )
    } else if (geom == "ridgeline") {
        p <- .plotReadsPerCellRidgeline(
            data = data,
            cutoff = cutoff,
            fill = fill
        )
    } else if (geom == "violin") {
        p <- .plotReadsPerCellViolin(
            data = data,
            cutoff = cutoff,
            fill = fill
        )
    }

    # Add title and subtitle containing cutoff information
    p <- p +
        labs(
            color = "",
            fill = "",
            title = "reads per cell",
            subtitle = paste("cutoff", cutoff, sep = " = ")
        )

    p
}



# Standard (raw) ---------------------------------------------------------------
.plotReadsPerCellECDF <- function(
    data,
    cutoff = 0L,
    color = scale_color_viridis(discrete = TRUE)
) {
    assert_is_a_number(cutoff)
    assertIsColorScaleDiscreteOrNULL(color)

    p <- ggplot(
        data = data,
        mapping = aes_string(
            x = "nCount",
            color = "interestingGroups"
        )
    ) +
        stat_ecdf(geom = "step") +
        scale_x_continuous(trans = "log10")

    # Cutoff lines
    if (cutoff > 0L) {
        p <- p + .qcCutoffLine(xintercept = cutoff)
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



.plotReadsPerCellViolin <- function(
    data,
    cutoff = 0L,
    fill = scale_fill_viridis(discrete = TRUE)
) {
    assert_is_a_number(cutoff)
    assertIsFillScaleDiscreteOrNULL(fill)

    p <- ggplot(
        data = data,
        mapping = aes_string(
            x = "sampleName",
            y = "nCount",
            fill = "interestingGroups"
        )
    ) +
        geom_violin(
            alpha = qcPlotAlpha,
            color = lineColor,
            scale = "count"
        ) +
        scale_y_continuous(trans = "log10") +
        .medianLabels(data, medianCol = "nCount", digits = 0L) +
        theme(axis.text.x = element_text(angle = 90L, hjust = 1L))

    # Cutoff lines
    if (cutoff > 0L) {
        p <- p + .qcCutoffLine(yintercept = cutoff)
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



.plotReadsPerCellRidgeline <- function(
    data,
    cutoff = 0L,
    fill = scale_fill_viridis(discrete = TRUE)
) {
    assert_is_a_number(cutoff)
    assertIsFillScaleDiscreteOrNULL(fill)

    p <- ggplot(
        data = data,
        mapping = aes_string(
            x = "nCount",
            y = "sampleName",
            fill = "interestingGroups"
        )
    ) +
        geom_density_ridges(
            alpha = qcPlotAlpha,
            color = lineColor,
            panel_scaling = TRUE,
            scale = 10L
        ) +
        scale_x_continuous(trans = "log10") +
        .medianLabels(data, medianCol = "nCount", digits = 0L)

    # Cutoff lines
    if (cutoff > 0L && length(cutoff)) {
        p <- p + .qcCutoffLine(xintercept = cutoff)
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



# Proportional -----------------------------------------------------------------
#' Proportional Cellular Barcodes Data
#'
#' Modified version of Allon Klein Lab MATLAB code.
#'
#' @author Rory Kirchner, Michael Steinbaugh
#' @keywords internal
#' @noRd
#'
#' @param data Cellular barcodes tibble containing the raw read counts.
#'
#' @return `grouped_df`, grouped by `sampleID`.
.proportionalReadsPerCell <- function(
    data,
    sampleData,
    breaks = 100L
) {
    assert_is_all_of(data, "grouped_df")
    assert_is_subset(c("nCount", "sampleID"), colnames(data))
    assert_is_integer(data[["nCount"]])
    assert_is_factor(data[["sampleID"]])
    assert_is_an_integer(breaks)
    mclapply(
        X = levels(data[["sampleID"]]),
        FUN = function(sampleID) {
            subset <- data[data[["sampleID"]] == sampleID, , drop = FALSE]
            # Histogram of log10-transformed counts
            h <- hist(
                x = log10(subset[["nCount"]]),
                n = breaks,
                plot = FALSE
            )
            # Klein Lab MATLAB code reference
            # counts: fLog
            # mids: xLog
            proportion <- h[["counts"]] * (10L ^ h[["mids"]]) /
                sum(h[["counts"]] * (10L ^ h[["mids"]]))
            tibble(
                "sampleID" = sampleID,
                "log10Count" = h[["mids"]],
                "proportion" = proportion
            )
        }
    ) %>%
        bind_rows() %>%
        mutate_if(is.character, as.factor) %>%
        group_by(!!sym("sampleID")) %>%
        left_join(sampleData, by = "sampleID")
}



.plotReadsPerCellHistogram <- function(
    data,
    cutoff = 0L,
    color = scale_color_viridis(discrete = TRUE)
) {
    assert_is_a_number(cutoff)
    assertIsColorScaleDiscreteOrNULL(color)

    p <- ggplot(
        data = data,
        mapping = aes_string(
            x = "log10Count",
            y = "proportion",
            color = "interestingGroups"
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
    if (cutoff > 0L) {
        p <- p + .qcCutoffLine(xintercept = log10(cutoff))
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
    .plotReadsPerCell
)



#' @rdname plotReadsPerCell
#' @export
setMethod(
    "plotReadsPerCell",
    signature("seurat"),
    .plotReadsPerCell
)
