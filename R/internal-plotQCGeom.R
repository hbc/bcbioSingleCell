# FIXME
# validMedianGeom <- c("violin", "boxplot", "ridgeline")
# validQCGeomFlip <- c("violin", "boxplot")


# FIXME Handle `flip = TRUE` for histogram

.plotQCGeom <- function(
    metrics,
    geom = c("violin", "boxplot", "histogram", "ridgeline"),
    ...
) {
    assert_is_data.frame(metrics)
    geom <- match.arg(geom)
    if (geom == "boxplot") {
        fxn <- .plotQCBoxplot
    } else if (geom == "histogram") {
        fxn <- .plotQCHistogram
    } else if (geom == "ridgeline") {
        fxn <- .plotQCRidgeline
    } else if (geom == "violin") {
        fxn <- .plotQCViolin
    }
    assert_is_function(fxn)
    fxn(metrics, ...)
}



# geom =========================================================================
.plotQCBoxplot <- function(
    metrics,
    metricCol,
    min = 0L,
    max = Inf,
    interestingGroups = "sampleName",
    flip = TRUE,
    fill = scale_fill_viridis(discrete = TRUE)
) {
    assert_is_data.frame(metrics)
    assert_is_a_string(metricCol)
    assert_all_are_non_negative(c(min, max))
    min <- min(min)
    max <- max(max)
    assertFormalInterestingGroups(metrics, interestingGroups)
    assert_is_a_bool(flip)
    assertIsFillScaleDiscreteOrNULL(fill)

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

    # Label interesting groups
    p <- p + labs(fill = paste(interestingGroups, collapse = ":\n"))

    # Color palette
    if (is(fill, "ScaleDiscrete")) {
        p <- p + fill
    }

    # Median labels
    p <- p + .medianLabels(metrics, medianCol = metricCol)

    # Facets
    facets <- NULL
    if (isTRUE(.checkAggregate(object))) {
        facets <- "sampleNameAggregate"
    }
    if (is.character(facets)) {
        p <- p + facet_wrap(facets = facets, scales = "free_y")
    }

    # Flip samples on y-axis
    if (isTRUE(flip)) {
        p <- p + coord_flip()
    }

    p
}



.plotQCHistogram <- function(
    metrics,
    metricCol,
    min = 0L,
    max = Inf,
    interestingGroups = "sampleName",
    fill = scale_fill_viridis(discrete = TRUE)
) {
    assert_is_data.frame(metrics)
    assert_is_a_string(metricCol)
    assert_all_are_non_negative(c(min, max))
    min <- min(min)
    max <- max(max)
    assertFormalInterestingGroups(metrics, interestingGroups)
    assertIsFillScaleDiscreteOrNULL(fill)

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

    # Label interesting groups
    p <- p + labs(fill = paste(interestingGroups, collapse = ":\n"))

    # Color palette
    if (!is.null(fill)) {
        p <- p + fill
    }

    # Facets
    facets <- NULL
    if (isTRUE(.checkAggregate(object))) {
        facets <- "sampleNameAggregate"
    }
    if (is.character(facets)) {
        p <- p + facet_wrap(facets = facets, scales = "free_y")
    }

    # Label interesting groups
    p <- p + labs(fill = paste(interestingGroups, collapse = ":\n"))

    # Color palette
    if (!is.null(fill)) {
        p <- p + fill
    }

    # Median labels
    p <- p + .medianLabels(metrics, medianCol = metricCol)

    # Facets
    facets <- NULL
    if (isTRUE(.checkAggregate(object))) {
        facets <- "sampleNameAggregate"
    }
    if (is.character(facets)) {
        p <- p + facet_wrap(facets = facets, scales = "free_y")
    }

    # Flip samples on y-axis
    if (isTRUE(flip)) {
        p <- p + coord_flip()
    }

    p
}



.plotQCRidgeline <- function(
    metrics,
    metricCol,
    min = 0L,
    max = Inf,
    interestingGroups = "sampleName",
    flip = TRUE,
    fill = scale_fill_viridis(discrete = TRUE)
) {
    assert_is_data.frame(metrics)
    assert_is_a_string(metricCol)
    assert_all_are_non_negative(c(min, max))
    min <- min(min)
    max <- max(max)
    assertFormalInterestingGroups(metrics, interestingGroups)
    assert_is_a_bool(flip)
    assertIsFillScaleDiscreteOrNULL(fill)

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

    # Label interesting groups
    p <- p + labs(fill = paste(interestingGroups, collapse = ":\n"))

    # Color palette
    if (!is.null(fill)) {
        p <- p + fill
    }

    # Median labels
    p <- p + .medianLabels(metrics, medianCol = metricCol)

    # Facets
    facets <- NULL
    if (isTRUE(.checkAggregate(object))) {
        facets <- "sampleNameAggregate"
    }
    if (is.character(facets)) {
        p <- p + facet_wrap(facets = facets, scales = "free_y")
    }

    # Flip samples on y-axis
    if (isTRUE(flip)) {
        p <- p + coord_flip()
    }

    p
}



.plotQCScatterplot <- function(
    metrics,
    xCol,
    yCol,
    interestingGroups = "sampleName",
    color = scale_color_viridis(discrete = TRUE)
) {
    assert_is_data.frame(metrics)
    assert_is_a_string(xCol)
    assert_is_a_string(yCol)

    p <- ggplot(
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

    # Label interesting groups
    p <- p + labs(color = paste(interestingGroups, collapse = ":\n"))

    # Color palette
    if (!is.null(fill)) {
        p <- p + fill
    }

    # Median labels
    p <- p + .medianLabels(metrics, medianCol = metricCol)

    # Facets
    facets <- NULL
    if (isTRUE(.checkAggregate(object))) {
        facets <- "sampleNameAggregate"
    }
    if (is.character(facets)) {
        p <- p + facet_wrap(facets = facets, scales = "free_y")
    }

    # Flip samples on y-axis
    if (isTRUE(flip)) {
        p <- p + coord_flip()
    }

    p
}



.plotQCViolin <- function(
    metrics,
    metricCol,
    min = 0L,
    max = Inf,
    interestingGroups = "sampleName",
    flip = TRUE,
    fill = scale_fill_viridis(discrete = TRUE)
) {
    assert_is_data.frame(metrics)
    assert_is_a_string(metricCol)
    assert_all_are_non_negative(c(min, max))
    min <- min(min)
    max <- max(max)
    assertFormalInterestingGroups(metrics, interestingGroups)
    assert_is_a_bool(flip)
    assertIsFillScaleDiscreteOrNULL(fill)

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

    # Label interesting groups
    p <- p + labs(fill = paste(interestingGroups, collapse = ":\n"))

    # Color palette
    if (!is.null(fill)) {
        p <- p + fill
    }

    # Median labels
    p <- p + .medianLabels(metrics, medianCol = metricCol)

    # Facets
    facets <- NULL
    if (isTRUE(.checkAggregate(object))) {
        facets <- "sampleNameAggregate"
    }
    if (is.character(facets)) {
        p <- p + facet_wrap(facets = facets, scales = "free_y")
    }

    # Flip samples on y-axis
    if (isTRUE(flip)) {
        p <- p + coord_flip()
    }

    p
}
