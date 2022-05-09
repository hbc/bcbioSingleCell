#' @name plotReadsPerCell
#' @author Michael Steinbaugh, Rory Kirchner
#' @inherit AcidGenerics::plotReadsPerCell
#' @note Updated 2021-09-10.
#'
#' @inheritParams AcidRoxygen::params
#' @param ... Additional arguments.
#'
#' @param cutoffLine `logical(1)`.
#' Include a line marking the cutoff.
#'
#' @examples
#' data(bcb)
#'
#' ## bcbioSingleCell ====
#' plotReadsPerCell(bcb, geom = "histogram")
#' plotReadsPerCell(bcb, geom = "ecdf")
NULL



#' Proportional cellular barcodes data
#'
#' Modified version of Allon Klein Lab MATLAB code.
#'
#' @author Michael Steinbaugh, Rory Kirchner
#' @keywords internal
#' @note Updated 2022-05-07.
#' @noRd
#'
#' @param data `DataFrame`.
#' Raw read counts per cellular barcode.
#' Return from `.rawMetrics()` function.
#'
#' @return `DataFrame`.
.proportionalReadsPerCell <-
    function(data,
             sampleData,
             breaks = 100L) {
        assert(
            requireNamespaces("graphics"),
            is(data, "DataFrame"),
            isSubset(c("nRead", "sampleId"), colnames(data)),
            is.integer(data[["nRead"]]),
            is.factor(data[["sampleId"]]),
            is(sampleData, "DataFrame"),
            isInt(breaks)
        )
        sampleData[["sampleId"]] <- as.factor(rownames(sampleData))
        samples <- levels(data[["sampleId"]])
        list <- DataFrameList(lapply(
            X = samples,
            FUN = function(sampleId) {
                keep <- which(data[["sampleId"]] == sampleId)
                subset <- data[keep, , drop = FALSE]
                ## Histogram of log10-transformed counts.
                h <- graphics::hist(
                    x = log10(subset[["nRead"]]),
                    n = breaks,
                    plot = FALSE
                )
                ## Klein Lab MATLAB code reference.
                ## counts: fLog; mids: xLog
                proportion <-
                    h[["counts"]] *
                        (10L^h[["mids"]]) /
                        sum(h[["counts"]] * (10L^h[["mids"]]))
                DataFrame(
                    "sampleId" = factor(sampleId),
                    "log10Read" = h[["mids"]],
                    "proportion" = proportion
                )
            }
        ))
        out <- unlist(list, recursive = FALSE, use.names = FALSE)
        out <- leftJoin(out, sampleData, by = "sampleId")
        out
    }



#' Plot proportional reads per cell histogram
#'
#' @note Updated 2022-05-07.
#' @noRd
#'
#' @param data Return from `.proportionalReadsPerCell()` function.
#'
#' @return `ggplot`.
.plotReadsPerCellHistogram <-
    function(data,
             min = 0L) {
        assert(is(data, "DataFrame"))
        p <- ggplot(
            data = as.data.frame(data),
            mapping = aes(
                x = !!sym("log10Read"),
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
        ## Cutoff line.
        if (min > 0L) {
            p <- p + acid_geom_abline(xintercept = log10(min))
        }
        ## Color palette.
        p <- p + autoDiscreteColorScale()
        ## Facets.
        facets <- NULL
        if (isSubset("aggregate", colnames(data))) {
            facets <- c(facets, "aggregate")
        }
        if (is.character(facets)) {
            p <- p + facet_wrap(facets = syms(facets), scales = "free")
        }
        ## Return.
        p
    }



## Updated 2022-05-07.
.plotReadsPerCellBoxplot <-
    function(data,
             min = 0L) {
        assert(is(data, "DataFrame"))
        p <- ggplot(
            data = as.data.frame(data),
            mapping = aes(
                x = !!sym("sampleName"),
                y = !!sym("nRead"),
                fill = !!sym("interestingGroups")
            )
        ) +
            geom_boxplot(color = "black", outlier.shape = NA) +
            scale_y_continuous(trans = "log10") +
            acid_geom_label_average(
                data = as.data.frame(data),
                col = "nRead",
                digits = 0L
            ) +
            labs(
                x = NULL,
                y = "reads per cell"
            )
        ## Cutoff line.
        if (min > 0L) {
            p <- p + acid_geom_abline(yintercept = min)
        }
        ## Color palette.
        p <- p + autoDiscreteFillScale()
        ## Facets.
        facets <- NULL
        if (isSubset("aggregate", colnames(data))) {
            facets <- c(facets, "aggregate")
        }
        if (is.character(facets)) {
            p <- p + facet_wrap(facets = syms(facets), scales = "free")
        }
        ## Return.
        p
    }



## Updated 2022-05-07.
.plotReadsPerCellECDF <-
    function(data,
             min = 0L) {
        assert(is(data, "DataFrame"))
        p <- ggplot(
            data = as.data.frame(data),
            mapping = aes(
                x = !!sym("nRead"),
                color = !!sym("interestingGroups")
            )
        ) +
            stat_ecdf(geom = "step", size = 1L) +
            labs(
                x = "reads per cell",
                y = "frequency"
            ) +
            scale_x_continuous(trans = "log10")
        ## Cutoff line.
        if (min > 0L) {
            p <- p + acid_geom_abline(xintercept = min)
        }
        ## Color palette.
        p <- p + autoDiscreteColorScale()
        ## Facets.
        facets <- NULL
        if (isSubset("aggregate", colnames(data))) {
            facets <- c(facets, "aggregate")
        }
        if (is.character(facets)) {
            p <- p + facet_wrap(facets = syms(facets), scales = "free")
        }
        ## Return.
        p
    }



## Updated 2022-05-07.
.plotReadsPerCellRidgeline <-
    function(data,
             min = 0L) {
        assert(is(data, "DataFrame"))
        p <- ggplot(
            data = as.data.frame(data),
            mapping = aes(
                x = !!sym("nRead"),
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
            acid_geom_label_average(
                data = as.data.frame(data),
                col = "nRead",
                digits = 0L
            ) +
            labs(
                x = "reads per cell",
                y = NULL
            )
        ## Cutoff line.
        if (min > 0L) {
            p <- p + acid_geom_abline(xintercept = min)
        }
        ## Color palette.
        p <- p + autoDiscreteFillScale()
        ## Facets.
        facets <- NULL
        if (isSubset("aggregate", colnames(data))) {
            facets <- c(facets, "aggregate")
        }
        if (is.character(facets)) {
            p <- p + facet_wrap(facets = syms(facets), scales = "free")
        }
        p
    }



## Updated 2022-05-07.
.plotReadsPerCellViolin <-
    function(data,
             min = 0L) {
        assert(is(data, "DataFrame"))
        p <- ggplot(
            data = as.data.frame(data),
            mapping = aes(
                x = !!sym("sampleName"),
                y = !!sym("nRead"),
                fill = !!sym("interestingGroups")
            )
        ) +
            geom_violin(
                color = "black",
                scale = "count"
            ) +
            scale_y_continuous(trans = "log10") +
            acid_geom_label_average(
                data = as.data.frame(data),
                col = "nRead",
                digits = 0L
            ) +
            labs(
                x = NULL,
                y = "reads per cell"
            )
        ## Cutoff line.
        if (min > 0L) {
            p <- p + acid_geom_abline(yintercept = min)
        }
        ## Color palette.
        p <- p + autoDiscreteFillScale()
        ## Facets.
        facets <- NULL
        if (isSubset("aggregate", colnames(data))) {
            facets <- c(facets, "aggregate")
        }
        if (is.character(facets)) {
            p <- p + facet_wrap(facets = syms(facets), scales = "free")
        }
        ## Return.
        p
    }



## Updated 2022-05-07.
`plotReadsPerCell,bcbioSingleCell` <- # nolint
    function(object,
             interestingGroups = NULL,
             geom,
             cutoffLine = FALSE,
             title = "Reads per cell") {
        validObject(object)
        assert(isString(title, nullOK = TRUE))
        interestingGroups(object) <-
            matchInterestingGroups(object, interestingGroups)
        interestingGroups <- interestingGroups(object)
        geom <- match.arg(geom)
        ## Minimum reads per barcode cutoff (for unfiltered data).
        if (!is.null(metadata(object)[["filterCells"]])) {
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
        ## This step will intentionally error for filtered objects.
        data <- .rawMetrics(object)
        p <- switch(
            EXPR = geom,
            boxplot = do.call(
                what = .plotReadsPerCellBoxplot,
                args = list(
                    "data" = data,
                    "min" = min
                )
            ),
            ecdf = do.call(
                what = .plotReadsPerCellECDF,
                args = list(
                    "data" = data,
                    "min" = min
                )
            ),
            histogram = {
                data <- do.call(
                    what = .proportionalReadsPerCell,
                    args = list(
                        "data" = data,
                        "sampleData" = sampleData(object)
                    )
                )
                do.call(
                    what = .plotReadsPerCellHistogram,
                    args = list(
                        "data" = data,
                        "min" = min
                    )
                )
            },
            ridgeline = do.call(
                what = .plotReadsPerCellRidgeline,
                args = list(
                    "data" = data,
                    "min" = min
                )
            ),
            violin = do.call(
                what = .plotReadsPerCellViolin,
                args = list(
                    "data" = data,
                    "min" = min
                )
            )
        )
        ## Add title and subtitle containing cutoff information.
        p <- p +
            labs(
                title = title,
                subtitle = subtitle,
                color = paste(interestingGroups, collapse = ":\n"),
                fill = paste(interestingGroups, collapse = ":\n")
            )
        ## Return.
        p
    }

formals(`plotReadsPerCell,bcbioSingleCell`)[["geom"]] <- .geom



#' @rdname plotReadsPerCell
#' @export
setMethod(
    f = "plotReadsPerCell",
    signature = signature(object = "bcbioSingleCell"),
    definition = `plotReadsPerCell,bcbioSingleCell`
)
