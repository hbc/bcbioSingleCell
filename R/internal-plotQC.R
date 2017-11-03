validQCGeom <- c("boxplot", "histogram", "ridgeline", "violin")
validGCGeomFlip <- c("boxplot", "violin")


.checkGeom <- function(geom) {
    if (!geom %in% validQCGeom) {
        stop(paste(
            "Valid formats:", toString(validQCGeom)
        ), call. = FALSE)
    }
}



#' @importFrom ggplot2 aes_string element_text geom_boxplot ggplot labs
#'   scale_y_sqrt theme
.plotQCBoxplot <- function(metrics, col, min, max) {
    p <- ggplot(
        metrics,
        mapping = aes_string(
            x = "sampleName",
            y = col,
            fill = "interestingGroups")
    ) +
        labs(x = "sample",
             y = "genes per cell") +
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
.plotQCHistogram <- function(metrics, col, min, max) {
    p <- ggplot(
        metrics,
        mapping = aes_string(
            x = col,
            fill = "interestingGroups")
    ) +
        labs(x = "genes per cell") +
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
.plotQCRidgeline <- function(metrics, col, min, max) {
    p <- ggplot(
        metrics,
        mapping = aes_string(
            x = col,
            y = "sampleName",
            fill = "interestingGroups")
    ) +
        labs(x = "genes per cell",
             y = "sample") +
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



#' @importFrom ggplot2 aes_string element_text geom_violin ggplot labs
#'   scale_y_sqrt theme
.plotQCViolin <- function(metrics, col, min, max) {
    p <- ggplot(
        metrics,
        mapping = aes_string(
            x = "sampleName",
            y = col,
            fill = "interestingGroups")
    ) +
        labs(x = "sample",
             y = "genes per cell") +
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
