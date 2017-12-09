#' Quality Control Plots
#'
#' Utility function that loops our standard quality control plots, for easy
#' visualization.
#'
#' @rdname plotQC
#' @name plotQC
#' @author Michael Steinbaugh
#'
#' @importFrom basejump plotQC
#'
#' @inheritParams AllGenerics
#' @inheritParams metrics
#'
#' @param geom Plot type. Supported formats: `boxplot`, `histogram`,
#'   `ridgeline`, and `violin`. Applies to [plotUMIsPerCell()],
#'   [plotGenesPerCell()], [plotMitoRatio()], and [plotNovelty()] output.
#' @param headerLevel R Markdown header level.
#' @param legend Include plot legend.
#' @param return
#'   - `grid`: [cowplot::plot_grid()] graphical output.
#'   - `list`: [ggplot] [list].
#'   - `markdown`: R Markdown report, with reports separated by headers.
#'
#' @return R Markdown template code for quality control analysis.
#'
#' @examples
#' load(system.file(
#'     file.path("extdata", "bcb.rda"),
#'     package = "bcbioSingleCell"))
#' load(system.file(
#'     file.path("extdata", "seurat.rda"),
#'     package = "bcbioSingleCell"))
#'
#' # bcbioSingleCell
#' plotQC(bcb)
#'
#' # seurat
#' plotQC(seurat)
NULL



# Constructors =================================================================
validMedianGeom <- c(
    "boxplot",
    "ridgeline",
    "violin")
validQCGeomFlip <- c(
    "boxplot",
    "violin")



#' @importFrom cowplot plot_grid
#' @importFrom ggplot2 theme
.plotQC <- function(
    object,
    interestingGroups,
    geom = "violin",
    headerLevel = 2,
    legend = FALSE,
    return = "grid") {
    if (missing(interestingGroups)) {
        interestingGroups <- basejump::interestingGroups(object)
    }
    plotlist <- list(
        plotReadsPerCell = plotReadsPerCell(
            object,
            interestingGroups = interestingGroups),
        plotCellCounts = plotCellCounts(
            object,
            interestingGroups = interestingGroups),
        plotUMIsPerCell = plotUMIsPerCell(
            object,
            interestingGroups = interestingGroups,
            geom = geom),
        plotGenesPerCell = plotGenesPerCell(
            object,
            interestingGroups = interestingGroups,
            geom = geom),
        plotUMIsVsGenes = plotUMIsVsGenes(
            object,
            interestingGroups = interestingGroups),
        plotMitoRatio = plotMitoRatio(
            object,
            interestingGroups = interestingGroups,
            geom = geom),
        plotNovelty = plotNovelty(
            object,
            interestingGroups = interestingGroups,
            geom = geom)
    )
    # Remove any `NULL`` plots. This is useful for nuking the
    # `plotReadsPerCell()` return on an object that doesn't contain raw cellular
    # barcode counts.
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
        plot_grid(plotlist = plotlist, labels = "AUTO")
    } else {
        mdHeader(
            "Quality control plots",
            level = headerLevel,
            tabset = TRUE,
            asis = TRUE)

        # Reads per cell currently only supported for bcbio runs
        if (!is.null(plotlist[["plotReadsPerCell"]])) {
            mdHeader(
                "Reads per cell",
                level = headerLevel + 1,
                asis = TRUE)
            show(plotlist[["plotReadsPerCell"]])
        }

        mdHeader(
            "Cell counts",
            level = headerLevel + 1,
            asis = TRUE)
        show(plotlist[["plotCellCounts"]])

        mdHeader(
            "UMI counts per cell",
            level = headerLevel + 1,
            asis = TRUE)
        show(plotlist[["plotUMIsPerCell"]])

        mdHeader(
            "Genes detected",
            level = headerLevel + 1,
            asis = TRUE)
        show(plotlist[["plotGenesPerCell"]])

        mdHeader(
            "UMIs vs. genes",
            level = headerLevel + 1,
            asis = TRUE)
        show(plotlist[["plotUMIsVsGenes"]])

        mdHeader(
            "Mitochondrial counts ratio",
            level = headerLevel + 1,
            asis = TRUE)
        show(plotlist[["plotMitoRatio"]])

        mdHeader(
            "Novelty",
            level = headerLevel + 1,
            asis = TRUE)
        show(plotlist[["plotNovelty"]])
    }
}



.plotQCGeom <- function(..., geom = "violin") {
    validGeom <- c(
        "boxplot",
        "histogram",
        "ridgeline",
        "violin")
    if (geom == "boxplot") {
        .plotQCBoxplot(...)
    } else if (geom == "histogram") {
        .plotQCHistogram(...)
    } else if (geom == "ridgeline") {
        .plotQCRidgeline(...)
    } else if (geom == "violin") {
        .plotQCViolin(...)
    } else {
        stop(paste(
            "Valid formats:", toString(validGeom)
        ), call. = FALSE)
    }
}



#' @importFrom ggplot2 aes_string element_text geom_boxplot ggplot labs
#'   scale_y_sqrt theme
.plotQCBoxplot <- function(metrics, metricCol, min = 0, max = Inf) {
    if (!is.numeric(min)) min <- 0
    if (!is.numeric(max)) max <- Inf
    p <- ggplot(
        metrics,
        mapping = aes_string(
            x = "sampleName",
            y = metricCol,
            fill = "interestingGroups")
    ) +
        geom_boxplot(color = lineColor, outlier.shape = NA) +
        scale_y_sqrt() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))

    # Cutoff lines
    if (min > 0) {
        p <- p + .qcCutoffLine(yintercept = min)
    }
    if (max < Inf) {
        p <- p + .qcCutoffLine(yintercept = max)
    }

    p
}



#' @importFrom ggplot2 aes_string element_text geom_histogram ggplot labs
#'   scale_x_sqrt scale_y_sqrt theme
.plotQCHistogram <- function(metrics, metricCol, min = 0, max = Inf) {
    if (!is.numeric(min)) min <- 0
    if (!is.numeric(max)) max <- Inf
    p <- ggplot(
        metrics,
        mapping = aes_string(
            x = metricCol,
            fill = "interestingGroups")
    ) +
        geom_histogram(bins = bins) +
        scale_x_sqrt() +
        scale_y_sqrt() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))

    # Cutoff lines
    if (min > 0) {
        p <- p + .qcCutoffLine(xintercept = min)
    }
    if (max < Inf) {
        p <- p + .qcCutoffLine(xintercept = max)
    }

    p
}



#' @importFrom ggplot2 aes_string element_text geom_boxplot ggplot labs
#'   scale_x_sqrt theme
#' @importFrom ggridges geom_density_ridges
.plotQCRidgeline <- function(metrics, metricCol, min = 0, max = Inf) {
    if (!is.numeric(min)) min <- 0
    if (!is.numeric(max)) max <- Inf
    p <- ggplot(
        metrics,
        mapping = aes_string(
            x = metricCol,
            y = "sampleName",
            fill = "interestingGroups")
    ) +
        geom_density_ridges(
            alpha = qcPlotAlpha,
            color = lineColor,
            panel_scaling = TRUE,
            scale = 10) +
        scale_x_sqrt() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))

    # Cutoff lines
    if (min > 0) {
        p <- p + .qcCutoffLine(xintercept = min)
    }
    if (max < Inf) {
        p <- p + .qcCutoffLine(xintercept = max)
    }

    p
}



#' @importFrom dplyr filter pull
#' @importFrom tibble rownames_to_column
#' @importFrom viridis scale_color_viridis
.plotQCScatterplot <- function(metrics, xCol, yCol) {
    ggplot(
        metrics,
        mapping = aes_string(
            x = xCol,
            y = yCol,
            color = "interestingGroups")
    ) +
        geom_point(alpha = 0.25, size = 0.8) +
        geom_smooth(se = FALSE, size = 1.5) +
        scale_x_sqrt() +
        scale_y_sqrt()
}



#' @importFrom ggplot2 aes_string element_text geom_violin ggplot labs
#'   scale_y_sqrt theme
.plotQCViolin <- function(metrics, metricCol, min = 0, max = Inf) {
    if (!is.numeric(min)) min <- 0
    if (!is.numeric(max)) max <- Inf
    p <- ggplot(
        metrics,
        mapping = aes_string(
            x = "sampleName",
            y = metricCol,
            fill = "interestingGroups")
    ) +
        geom_violin(
            color = lineColor,
            scale = "width",
            trim = TRUE) +
        scale_y_sqrt() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))

    # Cutoff lines
    if (min > 0) {
        p <- p + .qcCutoffLine(yintercept = min)
    }
    if (max < Inf) {
        p <- p + .qcCutoffLine(yintercept = max)
    }

    p
}



# Methods ======================================================================
#' @rdname plotQC
#' @export
setMethod(
    "plotQC",
    signature("bcbioSingleCell"),
    .plotQC)



#' @rdname plotQC
#' @export
setMethod(
    "plotQC",
    signature("seurat"),
    .plotQC)
