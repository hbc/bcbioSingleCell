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
#' plotReadsPerCell(indrops_small, geom = "histogram")
#' plotReadsPerCell(indrops_small, geom = "ecdf")
NULL



# Constructors =================================================================
.plotReadsPerCellBoxplot <- function(
    data,
    min = 0L,
    fill = getOption("bcbio.discrete.fill", NULL)
) {
    assert_is_data.frame(data)
    assertIsFillScaleDiscreteOrNULL(fill)

    p <- ggplot(
        data = data,
        mapping = aes(
            x = !!sym("sampleName"),
            y = !!sym("nCount"),
            fill = !!sym("interestingGroups")
        )
    ) +
        geom_boxplot(color = "black", outlier.shape = NA) +
        scale_y_continuous(trans = "log10") +
        bcbio_geom_label_average(data, col = "nCount", digits = 0L) +
        labs(
            x = NULL,
            y = "reads per cell"
        )

    # Cutoff line
    if (min > 0L) {
        p <- p + bcbio_geom_abline(yintercept = min)
    }

    # Color palette
    if (is(fill, "ScaleDiscrete")) {
        p <- p + fill
    }

    # Facets
    facets <- NULL
    if (.isAggregate(data)) {
        facets <- c(facets, "aggregate")
    }
    if (is.character(facets)) {
        p <- p + facet_wrap(facets = syms(facets), scales = "free")
    }

    p
}



.plotReadsPerCellECDF <- function(
    data,
    min = 0L,
    color = getOption("bcbio.discrete.color", NULL)
) {
    assert_is_data.frame(data)
    assertIsColorScaleDiscreteOrNULL(color)

    p <- ggplot(
        data = data,
        mapping = aes(
            x = !!sym("nCount"),
            color = !!sym("interestingGroups")
        )
    ) +
        stat_ecdf(geom = "step", size = 1L) +
        labs(
            x = "reads per cell",
            y = "frequency"
        ) +
        scale_x_continuous(trans = "log10")

    # Cutoff line
    if (min > 0L) {
        p <- p + bcbio_geom_abline(xintercept = min)
    }

    # Color palette
    if (is(color, "ScaleDiscrete")) {
        p <- p + color
    }

    # Facets
    facets <- NULL
    if (.isAggregate(data)) {
        facets <- c(facets, "aggregate")
    }
    if (is.character(facets)) {
        p <- p + facet_wrap(facets = syms(facets), scales = "free")
    }

    p
}



.plotReadsPerCellRidgeline <- function(
    data,
    min = 0L,
    fill = getOption("bcbio.discrete.fill", NULL)
) {
    assert_is_data.frame(data)
    assertIsFillScaleDiscreteOrNULL(fill)

    p <- ggplot(
        data = data,
        mapping = aes(
            x = !!sym("nCount"),
            y = !!sym("sampleName"),
            fill = !!sym("interestingGroups")
        )
    ) +
        geom_density_ridges(
            alpha = 0.75,
            color = "black",
            panel_scaling = TRUE,
            scale = 10L
        ) +
        scale_x_continuous(trans = "log10") +
        bcbio_geom_label_average(data, col = "nCount", digits = 0L) +
        labs(
            x = "reads per cell",
            y = NULL
        )

    # Cutoff line
    if (min > 0L) {
        p <- p + bcbio_geom_abline(xintercept = min)
    }

    # Color palette
    if (is(fill, "ScaleDiscrete")) {
        p <- p + fill
    }

    # Facets
    facets <- NULL
    if (.isAggregate(data)) {
        facets <- c(facets, "aggregate")
    }
    if (is.character(facets)) {
        p <- p + facet_wrap(facets = syms(facets), scales = "free")
    }

    p
}



.plotReadsPerCellViolin <- function(
    data,
    min = 0L,
    fill = getOption("bcbio.discrete.fill", NULL)
) {
    assert_is_data.frame(data)
    assertIsFillScaleDiscreteOrNULL(fill)

    p <- ggplot(
        data = data,
        mapping = aes(
            x = !!sym("sampleName"),
            y = !!sym("nCount"),
            fill = !!sym("interestingGroups")
        )
    ) +
        geom_violin(
            color = "black",
            scale = "count"
        ) +
        scale_y_continuous(trans = "log10") +
        bcbio_geom_label_average(data, col = "nCount", digits = 0L) +
        labs(
            x = NULL,
            y = "reads per cell"
        )

    # Cutoff line
    if (min > 0L) {
        p <- p + bcbio_geom_abline(yintercept = min)
    }

    # Color palette
    if (is(fill, "ScaleDiscrete")) {
        p <- p + fill
    }

    # Facets
    facets <- NULL
    if (.isAggregate(data)) {
        facets <- c(facets, "aggregate")
    }
    if (is.character(facets)) {
        p <- p + facet_wrap(facets = syms(facets), scales = "free")
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
    lapply(
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
                sampleID = sampleID,
                "log10Count" = h[["mids"]],
                proportion = proportion
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
    min = 0L,
    color = getOption("bcbio.discrete.color", NULL)
) {
    assert_is_data.frame(data)
    assertIsColorScaleDiscreteOrNULL(color)

    p <- ggplot(
        data = data,
        mapping = aes(
            x = !!sym("log10Count"),
            y = !!sym("proportion"),
            color = !!sym("interestingGroups")
        )
    ) +
        geom_step(
            alpha = 0.75,
            size = 1L
        ) +
        labs(
            x = "log10 reads per cell",
            y = "proportion of reads"
        )

    # Cutoff line
    if (min > 0L) {
        p <- p + bcbio_geom_abline(xintercept = log10(min))
    }

    # Color palette
    if (is(color, "ScaleDiscrete")) {
        p <- p + color
    }

    # Facets
    facets <- NULL
    if (.isAggregate(data)) {
        facets <- c(facets, "aggregate")
    }
    if (is.character(facets)) {
        p <- p + facet_wrap(facets = syms(facets), scales = "free")
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
        interestingGroups,
        geom = c("histogram", "ecdf", "violin", "ridgeline", "boxplot"),
        color = getOption("bcbio.discrete.color", NULL),
        fill = getOption("bcbio.discrete.fill", NULL),
        title = "reads per cell"
    ) {
        # Passthrough: color, fill
        validObject(object)
        interestingGroups <- .prepareInterestingGroups(
            object = object,
            interestingGroups = interestingGroups
        )
        geom <- match.arg(geom)
        assertIsAStringOrNULL(title)

        # Minimum reads per barcode cutoff (for unfiltered data)
        if (length(metadata(object)[["filterCells"]])) {
            min <- 0L
            subtitle <- NULL
        } else {
            min <- metadata(object)[["cellularBarcodeCutoff"]]
            subtitle <- paste("cutoff", min, sep = " = ")
        }
        assert_is_an_integer(min)

        # Obtain the sample metadata
        sampleData <- sampleData(
            object = object,
            interestingGroups = interestingGroups
        )
        assert_is_non_empty(sampleData)
        sampleData <- as.data.frame(sampleData)
        sampleData[["sampleID"]] <- as.factor(rownames(sampleData))

        # Obtain the read counts. Use the unfiltered reads stashed in the
        # metadata if available, otherwise use the metrics return.
        cbList <- metadata(object)[["cellularBarcodes"]]
        if (length(cbList)) {
            data <- .bindCellularBarcodes(cbList)
        } else {
            data <- metrics(object)
        }

        # Early return NULL if `nCount` isn't present
        if (!"nCount" %in% colnames(data)) {
            warning("object does not contain nCount column in `metrics()`")
            return(invisible())
        }

        data <- left_join(
            x = data[, c("sampleID", "nCount")],
            y = sampleData,
            by = "sampleID"
        ) %>%
            as_tibble() %>%
            group_by(!!sym("sampleID"))

        if (geom == "boxplot") {
            p <- .plotReadsPerCellBoxplot(
                data = data,
                fill = fill,
                min = min
            )
        } else if (geom == "ecdf") {
            p <- .plotReadsPerCellECDF(
                data = data,
                color = color,
                min = min
            )
        } else if (geom == "histogram") {
            sampleData <- sampleData(
                object = object,
                interestingGroups = interestingGroups
            ) %>%
                as.data.frame()
            sampleData[["sampleID"]] <- as.factor(rownames(sampleData))
            data <- .proportionalReadsPerCell(
                data = data,
                sampleData = sampleData
            )
            p <- .plotReadsPerCellHistogram(
                data = data,
                color = color,
                min = min
            )
        } else if (geom == "ridgeline") {
            p <- .plotReadsPerCellRidgeline(
                data = data,
                fill = fill,
                min = min
            )
        } else if (geom == "violin") {
            p <- .plotReadsPerCellViolin(
                data = data,
                fill = fill,
                min = min
            )
        }

        # Add title and subtitle containing cutoff information
        p <- p +
            labs(
                title = title,
                subtitle = subtitle,
                color = paste(interestingGroups, collapse = ":\n"),
                fill = paste(interestingGroups, collapse = ":\n")
            )

        p
    }
)
