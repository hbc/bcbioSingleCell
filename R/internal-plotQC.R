validQCGeom <- c("boxplot", "histogram", "ridgeline", "violin")
validGCGeomFlip <- c("boxplot", "violin")



.dynamicQCPlot <- function(..., geom) {
    if (!geom %in% validQCGeom) {
        stop(paste(
            "Valid formats:", toString(validQCGeom)
        ), call. = FALSE)
    }
    if (geom == "boxplot") {
        plotQCBoxplot(...)
    } else if (geom == "histogram") {
        .plotQCHistogram(...)
    } else if (geom == "ridgeline") {
        .plotQCRidgeline(...)
    } else if (geom == "violin") {
        .plotQCViolin(...)
    }
}



#' @importFrom ggplot2 aes_string element_text geom_boxplot ggplot labs
#'   scale_y_sqrt theme
.plotQCBoxplot <- function(metrics, col1, min, max) {
    p <- ggplot(
        metrics,
        mapping = aes_string(
            x = "sampleName",
            y = col1,
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
.plotQCHistogram <- function(metrics, col1, min, max) {
    p <- ggplot(
        metrics,
        mapping = aes_string(
            x = col1,
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
.plotQCRidgeline <- function(metrics, col1, min, max) {
    p <- ggplot(
        metrics,
        mapping = aes_string(
            x = col1,
            y = "sampleName",
            fill = "interestingGroups")
    ) +
        geom_density_ridges(
            alpha = qcPlotAlpha,
            color = lineColor,
            panel_scaling = TRUE,
            scale = qcRidgeScale) +
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
.plotQCScatterplot <- function(metrics, col1, col2) {
    ggplot(
        metrics,
        mapping = aes_string(
            x = col1,
            y = col2,
            color = "interestingGroups")
    ) +
        geom_point(
            alpha = 0.25,
            size = 0.8) +
        geom_smooth(
            method = "gam",
            se = FALSE,
            size = 1.5) +
        scale_x_sqrt() +
        scale_y_sqrt()
}



#' @importFrom ggplot2 aes_string element_text geom_violin ggplot labs
#'   scale_y_sqrt theme
.plotQCViolin <- function(metrics, col1, min, max) {
    p <- ggplot(
        metrics,
        mapping = aes_string(
            x = "sampleName",
            y = col1,
            fill = "interestingGroups")
    ) +
        geom_violin(
            alpha = qcPlotAlpha,
            color = lineColor,
            scale = "width") +
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
