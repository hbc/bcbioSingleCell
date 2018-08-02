#' Quality Control Plots
#'
#' Utility function that loops our standard quality control plots, for easy
#' visualization.
#'
#' @name plotQC
#' @family Quality Control Functions
#' @author Michael Steinbaugh
#'
#' @importFrom bcbioBase plotQC
#' @export
#'
#' @inheritParams general
#'
#' @return
#'   - `grid`: [cowplot::plot_grid()] graphical output.
#'   - `list`: `list` containing `ggplot` objects.
#'   - `markdown`: R Markdown report, with reports separated by headers.
#'
#' @examples
#' # bcbioSingleCell ====
#' plotQC(indrops_small)
#'
#' # SingleCellExperiment ====
#' plotQC(cellranger_small)
NULL



# Constructors =================================================================
# Plot a single quality control metric
.plotQCMetric <- function(
    object,
    metricCol,
    geom = c("violin", "ridgeline", "ecdf", "histogram", "boxplot"),
    interestingGroups,
    min = 0L,
    max = Inf,
    trans = "identity",
    ratio = FALSE,
    color = getOption("bcbio.discrete.color", NULL),
    fill = getOption("bcbio.discrete.fill", NULL),
    title = NULL
) {
    assert_is_a_string(metricCol)
    geom <- match.arg(geom)
    if (missing(interestingGroups)) {
        interestingGroups <- basejump::interestingGroups(object)
    } else {
        interestingGroups(object) <- interestingGroups
    }
    assert_is_character(interestingGroups)
    assert_all_are_non_negative(c(min, max))
    # Support for per sample filtering cutoffs
    min <- min(min)
    max <- max(max)
    if (isTRUE(ratio)) {
        assert_all_are_in_range(c(min, max), lower = 0L, upper = 1L)
    }
    assert_is_a_string(trans)
    assertIsFillScaleDiscreteOrNULL(fill)
    assertIsAStringOrNULL(title)

    data <- metrics(object, interestingGroups = interestingGroups)
    if (!metricCol %in% colnames(data)) {
        return(invisible())
    }

    mapping <- aes(
        color = !!sym("interestingGroups"),
        fill = !!sym("interestingGroups")
    )

    if (geom %in% c("boxplot", "violin")) {
        mapping[["x"]] <- as.symbol("sampleName")
        mapping[["y"]] <- as.symbol(metricCol)
    } else if (geom == "ridgeline") {
        # ridgeline flips the axes
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
            labs(x = NULL)
    } else if (geom == "ecdf") {
        p <- p +
            stat_ecdf(geom = "step", size = 1L) +
            scale_x_continuous(trans = trans) +
            labs(
                y = "frequency"
            )
    } else if (geom == "histogram") {
        p <- p +
            geom_histogram(
                bins = 200L,
                color = FALSE
            ) +
            scale_x_continuous(trans = trans) +
            scale_y_continuous()
    } else if (geom == "ridgeline") {
        p <- p +
            geom_density_ridges(
                alpha = 0.75,
                color = "black",
                panel_scaling = TRUE,
                scale = 10L
            ) +
            scale_x_continuous(trans = trans) +
            labs(y = NULL)
    } else if (geom == "violin") {
        p <- p +
            geom_violin(
                color = "black",
                scale = "area",
                trim = TRUE
            ) +
            scale_y_continuous(trans = trans) +
            labs(x = NULL)
    }

    # Cutoff lines
    if (geom %in% c("boxplot", "violin")) {
        if (min > 0L) {
            p <- p + bcbio_geom_abline(yintercept = min)
        }
        if (
            (max < Inf && identical(ratio, FALSE)) ||
            (max < 1L && identical(ratio, TRUE))
        ) {
            p <- p + bcbio_geom_abline(yintercept = max)
        }
    } else {
        if (min > 0L) {
            p <- p + bcbio_geom_abline(xintercept = min)
        }
        if (
            (max < Inf && identical(ratio, FALSE)) ||
            (max < 1L && identical(ratio, TRUE))
        ) {
            p <- p + bcbio_geom_abline(xintercept = max)
        }
    }

    # Label interesting groups
    p <- p +
        labs(
            title = title,
            color = paste(interestingGroups, collapse = ":\n"),
            fill = paste(interestingGroups, collapse = ":\n")
        )

    # Color palette
    if (geom == "ecdf") {
        if (is(color, "ScaleDiscrete")) {
            p <- p + color
        }
    } else {
        if (is(fill, "ScaleDiscrete")) {
            p <- p + fill
        }
    }

    # Median labels
    if (!geom %in% c("ecdf", "histogram")) {
        if (metricCol %in% c("log10GenesPerUMI", "mitoRatio")) {
            digits <- 2L
        } else {
            digits <- 0L
        }
        p <- p +
            bcbio_geom_label_average(data, col = metricCol, digits = digits)
    }

    # Facets
    facets <- NULL
    if (.isAggregate(data)) {
        facets <- "aggregate"
    }
    if (is.character(facets)) {
        p <- p + facet_wrap(facets = syms(facets), scales = "free")
    }

    p
}



# Compare two quality control metrics
.plotQCScatterplot <- function(
    object,
    xCol,
    yCol,
    xTrans = "identity",
    yTrans = "identity",
    interestingGroups,
    trendline = FALSE,
    color = getOption("bcbio.discrete.color", NULL),
    title = NULL
) {
    assert_is_a_string(xCol)
    assert_is_a_string(yCol)
    assert_is_a_string(xTrans)
    assert_is_a_string(yTrans)
    if (missing(interestingGroups)) {
        interestingGroups <- basejump::interestingGroups(object)
    } else {
        interestingGroups(object) <- interestingGroups
    }
    assert_is_character(interestingGroups)
    assertIsColorScaleDiscreteOrNULL(color)
    assertIsAStringOrNULL(title)

    data <- metrics(object, interestingGroups = interestingGroups)
    if (!all(c(xCol, yCol) %in% colnames(data))) {
        warning(paste(
            deparse(substitute(object)), "must contain",
            toString(c(xCol, yCol)),
            "columns in `metrics()`"
        ))
        return(invisible())
    }

    p <- ggplot(
        data = data,
        mapping = aes(
            x = !!sym(xCol),
            y = !!sym(yCol),
            color = !!sym("interestingGroups")
        )
    ) +
        geom_point(alpha = 0.5, size = 1L) +
        scale_x_continuous(trans = xTrans) +
        scale_y_continuous(trans = yTrans)

    if (isTRUE(trendline)) {
        # If `method = "gam"`, `mgcv` package is required.
        # Otherwise build checks will error.
        p <- p + geom_smooth(method = "glm", se = FALSE, size = 1L)
    }

    # Label interesting groups
    p <- p + labs(
        title = title,
        color = paste(interestingGroups, collapse = ":\n")
    )

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
#' @rdname plotQC
#' @export
setMethod(
    "plotQC",
    signature("SingleCellExperiment"),
    function(
        object,
        interestingGroups,
        geom = c("violin", "ridgeline", "ecdf", "histogram", "boxplot"),
        headerLevel = 2L,
        legend = getOption("bcbio.legend", FALSE),
        return = c("grid", "list", "markdown")
    ) {
        if (missing(interestingGroups)) {
            interestingGroups <- basejump::interestingGroups(object)
        } else {
            interestingGroups(object) <- interestingGroups
        }
        assert_is_character(interestingGroups)
        geom <- match.arg(geom)
        return <- match.arg(return)

        plotCellCounts <- plotCellCounts(object)
        plotReadsPerCell <- NULL
        plotUMIsPerCell <- plotUMIsPerCell(object, geom = geom)
        plotGenesPerCell <- plotGenesPerCell(object, geom = geom)
        plotUMIsVsGenes <- plotUMIsVsGenes(object)
        plotNovelty <- plotNovelty(object, geom = geom)
        plotMitoRatio <- plotMitoRatio(object, geom = geom)
        plotZerosVsDepth <- plotZerosVsDepth(object)

        if (is(object, "bcbioSingleCell")) {
            # Don't show cell counts for unfiltered bcbio datasets
            if (!length(metadata(object)[["filterCells"]])) {
                plotCellCounts <- NULL
            }
            # Raw read counts are only stashed in bcbioSingleCell objects
            plotReadsPerCell <- plotReadsPerCell(object, geom = geom)
        }

        plotlist <- list(
            "Cell Counts" = plotCellCounts,
            "Reads per Cell" = plotReadsPerCell,
            "UMIs per Cell" = plotUMIsPerCell,
            "Genes per Cell" = plotGenesPerCell,
            "UMIs vs. Genes" = plotUMIsVsGenes,
            "Novelty" = plotNovelty,
            "Mito Ratio" = plotMitoRatio,
            "Zeros vs. Depth" = plotZerosVsDepth
        )

        # Remove any `NULL` plots. This is useful for nuking the
        # `plotReadsPerCell()` return on an object that doesn't contain raw
        # cellular barcode counts.
        plotlist <- Filter(Negate(is.null), plotlist)

        # Hide the legends, if desired
        if (identical(legend, FALSE)) {
            .hideLegend <- function(gg) {
                gg + theme(legend.position = "none")
            }
            plotlist <- lapply(plotlist, .hideLegend)
        }

        # Grid return mode
        if (return == "list") {
            plotlist
        } else if (return == "grid") {
            plot_grid(plotlist = plotlist)
        } else if (return == "markdown") {
            markdownHeader(
                text = "Filtered quality control metrics",
                level = headerLevel,
                tabset = TRUE,
                asis = TRUE
            )
            markdownPlotlist(
                plotlist = plotlist,
                headerLevel = headerLevel + 1L
            )
        }
    }
)
