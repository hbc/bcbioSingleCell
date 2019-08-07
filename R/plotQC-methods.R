#' @name plotQC
#' @author Michael Steinbaugh
#' @include globals.R
#' @inherit acidplots::plotQC
#' @note Updated 2019-07-29.
#'
#' @inheritParams acidroxygen::params
#' @param headerLevel `integer(1)` (`1`-`7`).
#'   R Markdown header level.
#' @param legend `logical(1)`.
#'   Include plot legend.
#' @param ... Additional arguments.
#'
#' @return
#'   - `grid`: Graphical output.
#'   - `list`: `list` containing `ggplot` objects.
#'   - `markdown`: R Markdown report, with reports separated by headers.
#'
#' @examples
#' data(indrops)
#' plotQC(indrops)
NULL



#' @rdname plotQC
#' @name plotQC
#' @importFrom bioverbs plotQC
#' @usage plotQC(object, ...)
#' @export
NULL



## Plot a single quality control metric.
## Updated 2019-07-29.
.plotQCMetric <- function(
    object,
    metricCol,
    geom,
    interestingGroups = NULL,
    min = 0L,
    max = Inf,
    trans = "identity",
    ratio = FALSE,
    color,
    fill,
    title = NULL
) {
    validObject(object)
    assert(
        is(object, "SingleCellExperiment"),
        isString(metricCol),
        all(isNonNegative(c(min, max))),
        isString(trans),
        isGGScale(fill, scale = "discrete", aes = "fill", nullOK = TRUE),
        isString(title, nullOK = TRUE)
    )
    geom <- match.arg(geom)
    interestingGroups(object) <-
        matchInterestingGroups(object, interestingGroups)

    ## Support for per sample filtering cutoffs.
    min <- min(min)
    max <- max(max)
    if (isTRUE(ratio)) {
        assert(all(isInRange(c(min, max), lower = 0L, upper = 1L)))
    }

    data <- metrics(object)
    if (!isSubset(metricCol, colnames(data))) {
        stop(sprintf("`%s` is not defined in `colData()`.", metricCol))
    } else if (anyNA(data[[metricCol]])) {
        stop(sprintf("`%s` in `colData()` contains NA values.", metricCol))
    } else if (all(data[[metricCol]] == 0L)) {
        stop(sprintf("`%s` in `colData()` contains only zeros.", metricCol))
    }

    mapping <- aes(
        color = !!sym("interestingGroups"),
        fill = !!sym("interestingGroups")
    )

    if (geom %in% c("boxplot", "violin")) {
        mapping[["x"]] <- as.symbol("sampleName")
        mapping[["y"]] <- as.symbol(metricCol)
    } else if (geom == "ridgeline") {
        ## Ridgeline flips the axes.
        mapping[["x"]] <- as.symbol(metricCol)
        mapping[["y"]] <- as.symbol("sampleName")
    } else if (geom %in% c("ecdf", "histogram")) {
        mapping[["x"]] <- as.symbol(metricCol)
    }

    p <- ggplot(data = data, mapping = mapping)

    if (geom == "boxplot") {
        p <- p +
            geom_boxplot(color = "black", outlier.shape = NA) +
            scale_y_continuous(trans = trans) +
            labs(
                x = NULL,
                y = makeLabel(metricCol)
            )
    } else if (geom == "ecdf") {
        p <- p +
            stat_ecdf(geom = "step", size = 1L) +
            scale_x_continuous(trans = trans) +
            labs(
                x = makeLabel(metricCol),
                y = "frequency"
            )
    } else if (geom == "histogram") {
        p <- p +
            geom_histogram(
                bins = 200L,
                color = FALSE
            ) +
            scale_x_continuous(trans = trans) +
            scale_y_continuous() +
            labs(
                x = makeLabel(metricCol)
            )
    } else if (geom == "ridgeline") {
        p <- p +
            geom_density_ridges(
                alpha = 0.75,
                color = "black",
                panel_scaling = TRUE,
                scale = 10L
            ) +
            scale_x_continuous(trans = trans) +
            labs(
                x = makeLabel(metricCol),
                y = NULL
            )
    } else if (geom == "violin") {
        p <- p +
            geom_violin(
                color = "black",
                scale = "area",
                trim = TRUE
            ) +
            scale_y_continuous(trans = trans) +
            labs(
                x = NULL,
                y = makeLabel(metricCol)
            )
    }

    ## Cutoff lines.
    if (geom %in% c("boxplot", "violin")) {
        if (min > 0L) {
            p <- p + acid_geom_abline(yintercept = min)
        }
        if (
            (max < Inf && identical(ratio, FALSE)) ||
            (max < 1L && identical(ratio, TRUE))
        ) {
            p <- p + acid_geom_abline(yintercept = max)
        }
    } else {
        if (min > 0L) {
            p <- p + acid_geom_abline(xintercept = min)
        }
        if (
            (max < Inf && identical(ratio, FALSE)) ||
            (max < 1L && identical(ratio, TRUE))
        ) {
            p <- p + acid_geom_abline(xintercept = max)
        }
    }

    ## Label interesting groups.
    p <- p +
        labs(
            title = title,
            color = paste(interestingGroups, collapse = ":\n"),
            fill = paste(interestingGroups, collapse = ":\n")
        )

    ## Color palette.
    if (geom == "ecdf") {
        if (is(color, "ScaleDiscrete")) {
            p <- p + color
        }
    } else {
        if (is(fill, "ScaleDiscrete")) {
            p <- p + fill
        }
    }

    ## Median labels.
    if (!geom %in% c("ecdf", "histogram")) {
        if (metricCol %in% c("log10GenesPerUMI", "mitoRatio")) {
            digits <- 2L
        } else {
            digits <- 0L
        }
        p <- p +
            acid_geom_label_average(data, col = metricCol, digits = digits)
    }

    ## Facets.
    facets <- NULL
    if (.isAggregate(data)) {
        facets <- "aggregate"
    }
    if (is.character(facets)) {
        p <- p + facet_wrap(facets = syms(facets), scales = "free")
    }

    p
}

formals(`.plotQCMetric`)[c("color", "fill")] <-
    formalsList[c("color.discrete", "fill.discrete")]
formals(`.plotQCMetric`)[["geom"]] <- geom



## Compare two quality control metrics.
## Updated 2019-07-27.
.plotQCScatterplot <- function(
    object,
    xCol,
    yCol,
    xTrans = "identity",
    yTrans = "identity",
    interestingGroups = NULL,
    trendline = FALSE,
    color = getOption("basejump.discrete.color", NULL),
    title = NULL
) {
    validObject(object)
    assert(
        is(object, "SingleCellExperiment"),
        isString(xCol),
        isString(yCol),
        isString(xTrans),
        isString(yTrans),
        isGGScale(color, scale = "discrete", aes = "colour", nullOK = TRUE),
        isString(title, nullOK = TRUE)
    )
    interestingGroups(object) <-
        matchInterestingGroups(object, interestingGroups)

    data <- metrics(object)
    if (!isSubset(c(xCol, yCol), colnames(data))) {
        stop(sprintf(
            "%s are not defined in `colData()`.",
            toString(c(xCol, yCol))
        ))
    } else if (anyNA(data[[xCol]])) {
        stop(sprintf("`%s` in `colData()` contains NA values.", xCol))
    } else if (anyNA(data[[yCol]])) {
        stop(sprintf("`%s` in `colData()` contains NA values.", yCol))
    } else if (all(data[[xCol]] == 0L)) {
        stop(sprintf("`%s` in `colData()` contains only zeros.", xCol))
    } else if (all(data[[yCol]] == 0L)) {
        stop(sprintf("`%s` in `colData()` contains only zeros.", yCol))
    }

    p <- ggplot(
        data = data,
        mapping = aes(
            x = !!sym(xCol),
            y = !!sym(yCol),
            colour = !!sym("interestingGroups")
        )
    ) +
        geom_point(alpha = 0.5, size = 1L) +
        scale_x_continuous(trans = xTrans) +
        scale_y_continuous(trans = yTrans)

    if (isTRUE(trendline)) {
        ## If `method = "gam"`, `mgcv` package is required.
        ## Otherwise build checks will error.
        p <- p + geom_smooth(method = "glm", se = FALSE, size = 1L)
    }

    ## Set the labels.
    p <- p + labs(
        x = makeLabel(xCol),
        y = makeLabel(yCol),
        title = title,
        colour = paste(interestingGroups, collapse = ":\n")
    )

    ## Color palette.
    if (is(color, "ScaleDiscrete")) {
        p <- p + color
    }

    ## Facets.
    facets <- NULL
    if (.isAggregate(data)) {
        facets <- c(facets, "aggregate")
    }
    if (is.character(facets)) {
        p <- p + facet_wrap(facets = syms(facets), scales = "free")
    }

    p
}



## Updated 2019-07-24.
`plotQC,bcbioSingleCell` <-  # nolint
    function(
        object,
        interestingGroups = NULL,
        geom,
        headerLevel = 2L,
        legend,
        return = c("grid", "list", "markdown")
    ) {
        validObject(object)
        assert(isHeaderLevel(headerLevel))
        interestingGroups(object) <-
            matchInterestingGroups(object, interestingGroups)
        geom <- match.arg(geom)
        return <- match.arg(return)

        ## Don't show cell counts for unfiltered datasets.
        if (!is.null(metadata(object)[["filterCells"]])) {
            plotCellCounts <- plotCellCounts(object)
            plotZerosVsDepth <- NULL
        } else {
            plotCellCounts <- NULL
            plotZerosVsDepth <- plotZerosVsDepth(object)
        }

        plotUMIsPerCell <- plotUMIsPerCell(object, geom = geom)
        plotGenesPerCell <- plotGenesPerCell(object, geom = geom)
        plotUMIsVsGenes <- plotUMIsVsGenes(object)
        plotNovelty <- plotNovelty(object, geom = geom)
        plotMitoRatio <- plotMitoRatio(object, geom = geom)

        list <- list(
            "Cell counts" = plotCellCounts,
            "UMIs per cell" = plotUMIsPerCell,
            "Genes per cell" = plotGenesPerCell,
            "UMIs vs. genes" = plotUMIsVsGenes,
            "Novelty" = plotNovelty,
            "Mito ratio" = plotMitoRatio,
            "Zeros vs. depth" = plotZerosVsDepth
        )

        ## Remove any `NULL` plots. This is useful for nuking the
        ## `plotReadsPerCell` return on an object that doesn't contain raw
        ## cellular barcode counts.
        list <- Filter(f = Negate(is.null), x = list)

        ## Consistently show n plots.
        n <- 6L
        assert(hasLength(list, n = n))

        ## Hide the legends, if desired.
        if (identical(legend, FALSE)) {
            .hideLegend <- function(gg) {
                gg + theme(legend.position = "none")
            }
            list <- lapply(list, .hideLegend)
        }

        ## Return.
        if (return == "list") {
            names(list) <- camelCase(names(list))
            list
        } else if (return == "grid") {
            plot_grid(
                plotlist = list,
                ncol = n / 2L,
                nrow = 2L
            )
        } else if (return == "markdown") {
            markdownHeader(
                text = "Quality control metrics",
                level = headerLevel,
                tabset = TRUE,
                asis = TRUE
            )
            markdownPlots(list = list, headerLevel = headerLevel + 1L)
        }
    }

formals(`plotQC,bcbioSingleCell`)[["geom"]] <- geom
formals(`plotQC,bcbioSingleCell`)[["legend"]] <- formalsList[["legend"]]



#' @rdname plotQC
#' @export
setMethod(
    f = "plotQC",
    signature = signature("bcbioSingleCell"),
    definition = `plotQC,bcbioSingleCell`
)
