.plotQCGeom <- function(
    metrics,
    geom = c("violin", "boxplot", "histogram", "ridgeline"),
    ...
) {
    assert_is_data.frame(metrics)
    geom <- match.arg(geom)
    if (geom == "boxplot") {
        f <- .plotQCBoxplot
    } else if (geom == "histogram") {
        f <- .plotQCHistogram
    } else if (geom == "ridgeline") {
        f <- .plotQCRidgeline
    } else if (geom == "violin") {
        f <- .plotQCViolin
    }
    assert_is_function(f)
    f(metrics, ...)
}



# geom =========================================================================
.plotQCBoxplot <- function(
    metrics,
    metricCol,
    min = 0L,
    max = Inf
) {
    assert_is_data.frame(metrics)
    assert_is_a_string(metricCol)
    assert_all_are_non_negative(c(min, max))
    min <- min(min)
    max <- max(max)

    p <- ggplot(
        data = metrics,
        mapping = aes_string(
            x = "sampleName",
            y = metricCol,
            fill = "interestingGroups"
        )
    ) +
        geom_boxplot(color = lineColor, outlier.shape = NA) +
        scale_y_sqrt() +
        theme(axis.text.x = element_text(angle = 90L, hjust = 1L))

    # Cutoff lines
    if (min > 0L) {
        p <- p + .qcCutoffLine(yintercept = min)
    }
    if (max < Inf) {
        p <- p + .qcCutoffLine(yintercept = max)
    }

    p
}



.plotQCHistogram <- function(
    metrics,
    metricCol,
    min = 0L,
    max = Inf
) {
    assert_is_data.frame(metrics)
    assert_is_a_string(metricCol)
    assert_all_are_non_negative(c(min, max))
    min <- min(min)
    max <- max(max)

    p <- ggplot(
        data = metrics,
        mapping = aes_string(
            x = metricCol,
            fill = "interestingGroups"
        )
    ) +
        geom_histogram(bins = bins) +
        scale_x_sqrt() +
        scale_y_sqrt() +
        theme(axis.text.x = element_text(angle = 90L, hjust = 1L))

    # Cutoff lines
    if (min > 0L) {
        p <- p + .qcCutoffLine(xintercept = min)
    }
    if (max < Inf) {
        p <- p + .qcCutoffLine(xintercept = max)
    }

    p
}



.plotQCRidgeline <- function(
    metrics,
    metricCol,
    min = 0L,
    max = Inf
) {
    assert_is_data.frame(metrics)
    assert_is_a_string(metricCol)
    assert_all_are_non_negative(c(min, max))
    min <- min(min)
    max <- max(max)

    p <- ggplot(
        data = metrics,
        mapping = aes_string(
            x = metricCol,
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
        scale_x_sqrt() +
        theme(axis.text.x = element_text(angle = 90L, hjust = 1L))

    # Cutoff lines
    if (min > 0L) {
        p <- p + .qcCutoffLine(xintercept = min)
    }
    if (max < Inf) {
        p <- p + .qcCutoffLine(xintercept = max)
    }

    p
}



.plotQCScatterplot <- function(metrics, xCol, yCol) {
    assert_is_data.frame(metrics)
    assert_is_a_string(xCol)
    assert_is_a_string(yCol)

    ggplot(
        data = metrics,
        mapping = aes_string(
            x = xCol,
            y = yCol,
            color = "interestingGroups"
        )
    ) +
        geom_point(alpha = 0.25, size = 0.8) +
        # If `method = "gam"`, `mgcv` package is required.
        # Otherwise build checks will error.
        geom_smooth(
            method = "glm",
            se = FALSE,
            size = 1.5
        ) +
        scale_x_sqrt() +
        scale_y_sqrt()
}



.plotQCViolin <- function(
    metrics,
    metricCol,
    min = 0L,
    max = Inf
) {
    assert_is_data.frame(metrics)
    assert_is_a_string(metricCol)
    assert_all_are_non_negative(c(min, max))
    min <- min(min)
    max <- max(max)

    p <- ggplot(
        data = metrics,
        mapping = aes_string(
            x = "sampleName",
            y = metricCol,
            fill = "interestingGroups"
        )
    ) +
        geom_violin(
            color = lineColor,
            scale = "width",
            trim = TRUE
        ) +
        scale_y_sqrt() +
        theme(axis.text.x = element_text(angle = 90L, hjust = 1L))

    # Cutoff lines
    if (min > 0L) {
        p <- p + .qcCutoffLine(yintercept = min)
    }
    if (max < Inf) {
        p <- p + .qcCutoffLine(yintercept = max)
    }

    p
}
