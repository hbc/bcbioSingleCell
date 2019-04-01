#' @name plotReadsPerCell
#' @author Michael Steinbaugh, Rory Kirchner
#' @include globals.R
#' @inherit bioverbs::plotReadsPerCell
#' @inheritParams minimalism::params
#' @inheritParams basejump::params

#' @param cutoffLine `logical(1)`.
#'   Include a line marking the cutoff.
#'
#' @examples
#' data(indrops)
#' plotReadsPerCell(indrops, geom = "histogram")
#' plotReadsPerCell(indrops, geom = "ecdf")
NULL



#' @importFrom bioverbs plotReadsPerCell
#' @aliases NULL
#' @export
bioverbs::plotReadsPerCell



# Histogram ====================================================================
#' Proportional Cellular Barcodes Data
#'
#' Modified version of Allon Klein Lab MATLAB code.
#'
#' @author Rory Kirchner, Michael Steinbaugh
#' @keywords internal
#' @noRd
#'
#' @param data `tbl_df`. Raw read counts per cellular barcode.
#'
#' @return `grouped_df`. Grouped by `sampleID`.
.proportionalReadsPerCell <- function(
    data,
    sampleData,
    breaks = 100L
) {
    assert(
        is(data, "grouped_df"),
        isSubset(c("nCount", "sampleID"), colnames(data)),
        is.integer(data[["nCount"]]),
        is.factor(data[["sampleID"]]),
        is(sampleData, "DataFrame"),
        isInt(breaks)
    )

    sampleData <- sampleData %>%
        as_tibble(rownames = "sampleID") %>%
        mutate_all(as.factor)

    list <- lapply(
        X = levels(data[["sampleID"]]),
        FUN = function(sampleID) {
            subset <- data[data[["sampleID"]] == sampleID, , drop = FALSE]
            # Histogram of log10-transformed counts.
            h <- hist(
                x = log10(subset[["nCount"]]),
                n = breaks,
                plot = FALSE
            )
            # Klein Lab MATLAB code reference.
            # counts: fLog; mids: xLog
            proportion <- h[["counts"]] * (10L ^ h[["mids"]]) /
                sum(h[["counts"]] * (10L ^ h[["mids"]]))
            tibble(
                sampleID = sampleID,
                "log10Count" = h[["mids"]],
                proportion = proportion
            )
        }
    )

    list %>%
        bind_rows() %>%
        mutate_if(is.character, as.factor) %>%
        group_by(!!sym("sampleID")) %>%
        left_join(sampleData, by = "sampleID")
}



.plotReadsPerCellHistogram <- function(
    data,
    min = 0L,
    color = getOption("basejump.discrete.color", NULL)
) {
    assert(
        is.data.frame(data),
        isGGScale(color, scale = "discrete", aes = "colour", nullOK = TRUE)
    )

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

    # Cutoff line.
    if (min > 0L) {
        p <- p + acid_geom_abline(xintercept = log10(min))
    }

    # Color palette.
    if (is(color, "ScaleDiscrete")) {
        p <- p + color
    }

    # Facets.
    facets <- NULL
    if (.isAggregate(data)) {
        facets <- c(facets, "aggregate")
    }
    if (is.character(facets)) {
        p <- p + facet_wrap(facets = syms(facets), scales = "free")
    }

    p
}



# Boxplot ======================================================================
.plotReadsPerCellBoxplot <- function(
    data,
    min = 0L,
    fill = getOption("basejump.discrete.fill", NULL)
) {
    assert(
        is.data.frame(data),
        isGGScale(fill, scale = "discrete", aes = "fill", nullOK = TRUE)
    )

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
        acid_geom_label_average(data, col = "nCount", digits = 0L) +
        labs(
            x = NULL,
            y = "reads per cell"
        )

    # Cutoff line.
    if (min > 0L) {
        p <- p + acid_geom_abline(yintercept = min)
    }

    # Color palette.
    if (is(fill, "ScaleDiscrete")) {
        p <- p + fill
    }

    # Facets.
    facets <- NULL
    if (.isAggregate(data)) {
        facets <- c(facets, "aggregate")
    }
    if (is.character(facets)) {
        p <- p + facet_wrap(facets = syms(facets), scales = "free")
    }

    p
}



# ECDF =========================================================================
.plotReadsPerCellECDF <- function(
    data,
    min = 0L,
    color = getOption("basejump.discrete.color", NULL)
) {
    assert(
        is.data.frame(data),
        isGGScale(color, scale = "discrete", aes = "colour", nullOK = TRUE)
    )

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

    # Cutoff line.
    if (min > 0L) {
        p <- p + acid_geom_abline(xintercept = min)
    }

    # Color palette.
    if (is(color, "ScaleDiscrete")) {
        p <- p + color
    }

    # Facets.
    facets <- NULL
    if (.isAggregate(data)) {
        facets <- c(facets, "aggregate")
    }
    if (is.character(facets)) {
        p <- p + facet_wrap(facets = syms(facets), scales = "free")
    }

    p
}



# Ridgeline ====================================================================
.plotReadsPerCellRidgeline <- function(
    data,
    min = 0L,
    fill = getOption("basejump.discrete.fill", NULL)
) {
    assert(
        is.data.frame(data),
        isGGScale(fill, scale = "discrete", aes = "fill", nullOK = TRUE)
    )

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
        acid_geom_label_average(data, col = "nCount", digits = 0L) +
        labs(
            x = "reads per cell",
            y = NULL
        )

    # Cutoff line.
    if (min > 0L) {
        p <- p + acid_geom_abline(xintercept = min)
    }

    # Color palette.
    if (is(fill, "ScaleDiscrete")) {
        p <- p + fill
    }

    # Facets.
    facets <- NULL
    if (.isAggregate(data)) {
        facets <- c(facets, "aggregate")
    }
    if (is.character(facets)) {
        p <- p + facet_wrap(facets = syms(facets), scales = "free")
    }

    p
}



# Violin =======================================================================
.plotReadsPerCellViolin <- function(
    data,
    min = 0L,
    fill = getOption("basejump.discrete.fill", NULL)
) {
    assert(
        is.data.frame(data),
        isGGScale(fill, scale = "discrete", aes = "fill", nullOK = TRUE)
    )

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
        acid_geom_label_average(data, col = "nCount", digits = 0L) +
        labs(
            x = NULL,
            y = "reads per cell"
        )

    # Cutoff line.
    if (min > 0L) {
        p <- p + acid_geom_abline(yintercept = min)
    }

    # Color palette.
    if (is(fill, "ScaleDiscrete")) {
        p <- p + fill
    }

    # Facets.
    facets <- NULL
    if (.isAggregate(data)) {
        facets <- c(facets, "aggregate")
    }
    if (is.character(facets)) {
        p <- p + facet_wrap(facets = syms(facets), scales = "free")
    }

    p
}



# bcbioSingleCell ==============================================================
plotReadsPerCell.bcbioSingleCell <-  # nolint
    function(
        object,
        interestingGroups = NULL,
        geom,
        cutoffLine = FALSE,
        color,
        fill,
        title = "reads per cell"
    ) {
        # Passthrough: color, fill.
        validObject(object)
        assert(isString(title, nullOK = TRUE))
        interestingGroups(object) <-
            matchInterestingGroups(object, interestingGroups)
        geom <- match.arg(geom)

        # Minimum reads per barcode cutoff (for unfiltered data).
        if (length(metadata(object)[["filterCells"]])) {
            min <- 0L
            subtitle <- NULL
        } else {
            cutoff <- metadata(object)[["cellularBarcodeCutoff"]]
            subtitle <- paste("cutoff", cutoff, sep = " = ")
            if (isTRUE(cutoffLine)) {
                min <- cutoff
            } else {
                min <- 0L
            }
        }
        assert(isInt(min))

        data <- .rawMetrics(object)

        if (geom == "boxplot") {
            p <- do.call(
                what = .plotReadsPerCellBoxplot,
                args = list(
                    data = data,
                    fill = fill,
                    min = min
                )
            )
        } else if (geom == "ecdf") {
            p <- do.call(
                what = .plotReadsPerCellECDF,
                args = list(
                    data = data,
                    color = color,
                    min = min
                )
            )
        } else if (geom == "histogram") {
            data <- do.call(
                what = .proportionalReadsPerCell,
                args = list(
                    data = data,
                    sampleData = sampleData(object)
                )
            )
            p <- do.call(
                what = .plotReadsPerCellHistogram,
                args = list(
                    data = data,
                    color = color,
                    min = min
                )
            )
        } else if (geom == "ridgeline") {
            p <- do.call(
                what = .plotReadsPerCellRidgeline,
                args = list(
                    data = data,
                    fill = fill,
                    min = min
                )
            )
        } else if (geom == "violin") {
            p <- do.call(
                what = .plotReadsPerCellViolin,
                args = list(
                    data = data,
                    fill = fill,
                    min = min
                )
            )
        }

        # Add title and subtitle containing cutoff information.
        p <- p +
            labs(
                title = title,
                subtitle = subtitle,
                color = paste(interestingGroups, collapse = ":\n"),
                fill = paste(interestingGroups, collapse = ":\n")
            )

        p
    }

formals(plotReadsPerCell.bcbioSingleCell)[["color"]] <-
    formalsList[["color.discrete"]]
formals(plotReadsPerCell.bcbioSingleCell)[["fill"]] <-
    formalsList[["fill.discrete"]]
formals(plotReadsPerCell.bcbioSingleCell)[["geom"]] <- geom



#' @rdname plotReadsPerCell
#' @export
setMethod(
    f = "plotReadsPerCell",
    signature = signature("bcbioSingleCell"),
    definition = plotReadsPerCell.bcbioSingleCell
)
