# Plot a single quality control metric
.plotQCMetric <- function(
    object,
    metricCol,
    geom = c("boxplot", "ecdf", "histogram", "ridgeline", "violin"),
    interestingGroups,
    min = 0L,
    max = Inf,
    trans = "identity",
    color = scale_color_viridis(discrete = TRUE),
    fill = scale_fill_viridis(discrete = TRUE)
) {
    assert_is_a_string(metricCol)
    geom <- match.arg(geom)
    if (missing(interestingGroups)) {
        interestingGroups <- bcbioBase::interestingGroups(object)
    }
    if (missing(min)) {
        min <- metadata(object)[["filterParams"]][["minGenes"]]
        if (!is.numeric(min)) {
            min <- 0L
        }
    }
    if (missing(max)) {
        max <- metadata(object)[["filterParams"]][["maxGenes"]]
        if (!is.numeric(max)) {
            max <- Inf
        }
    }
    assert_all_are_non_negative(c(min, max))
    # Support for per sample filtering cutoffs
    min <- min(min)
    max <- max(max)
    assert_is_a_string(trans)
    assertIsFillScaleDiscreteOrNULL(fill)

    metrics <- metrics(object, interestingGroups = interestingGroups)
    if (!metricCol %in% colnames(metrics)) {
        warning(paste(
            deparse(substitute(object)),
            "does not contain", metricCol, "column in `metrics()`"
        ))
        return(invisible())
    }

    mapping <- aes_string(
        color = "interestingGroups",
        fill = "interestingGroups"
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

    p <- ggplot(data = metrics, mapping = mapping)

    if (geom == "boxplot") {
        p <- p +
            geom_boxplot(color = lineColor, outlier.shape = NA) +
            scale_y_continuous(trans = trans)
    } else if (geom == "ecdf") {
        p <- p +
            stat_ecdf(geom = "step", size = 1L) +
            scale_x_continuous(trans = trans) +
            labs(y = "frequency")
    } else if (geom == "histogram") {
        p <- p +
            geom_histogram(bins = bins) +
            scale_x_continuous(trans = trans) +
            scale_y_continuous(trans = trans)
    } else if (geom == "ridgeline") {
        p <- p +
            geom_density_ridges(
                alpha = qcPlotAlpha,
                color = lineColor,
                panel_scaling = TRUE,
                scale = 10L
            ) +
            scale_x_continuous(trans = trans)
    } else if (geom == "violin") {
        p <- p +
            geom_violin(
                color = lineColor,
                scale = "count",
                trim = TRUE
            ) +
            scale_y_continuous(trans = trans)
    }

    # Rotate the x-axis 90 degrees
    p <- p + theme(axis.text.x = element_text(angle = 90L, hjust = 1L))

    # Cutoff lines
    if (geom %in% c("boxplot", "violin")) {
        if (min > 0L) {
            p <- p + .qcCutoffLine(yintercept = min)
        }
        if (max < Inf) {
            p <- p + .qcCutoffLine(yintercept = max)
        }
    } else if (geom %in% c("histogram", "ridgeline")) {
        if (min > 0L) {
            p <- p + .qcCutoffLine(xintercept = min)
        }
        if (max < Inf) {
            p <- p + .qcCutoffLine(xintercept = max)
        }
    }

    # Label interesting groups
    p <- p +
        labs(
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
            .medianLabels(
                metrics,
                medianCol = metricCol,
                digits = digits
            )
    }

    # Facets
    facets <- NULL
    if (.isAggregate(object)) {
        facets <- "sampleNameAggregate"
    }
    if (is.character(facets)) {
        p <- p + facet_wrap(facets = facets, scales = "free")
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
    color = scale_color_viridis(discrete = TRUE)
) {
    if (missing(interestingGroups)) {
        interestingGroups <- bcbioBase::interestingGroups(object)
    }
    assert_is_a_string(xCol)
    assert_is_a_string(yCol)
    assert_is_a_string(xTrans)
    assert_is_a_string(yTrans)
    assertIsColorScaleDiscreteOrNULL(color)

    metrics <- metrics(object, interestingGroups = interestingGroups)
    if (!all(c(xCol, yCol) %in% colnames(metrics))) {
        warning(paste(
            deparse(substitute(object)), "must contain",
            toString(c(xCol, yCol)),
            "columns in `metrics()`"
        ))
        return(invisible())
    }

    p <- ggplot(
        data = metrics,
        mapping = aes_string(
            x = xCol,
            y = yCol,
            color = "interestingGroups"
        )
    ) +
        geom_point(alpha = 0.5, size = 1L) +
        # If `method = "gam"`, `mgcv` package is required.
        # Otherwise build checks will error.
        geom_smooth(method = "glm", se = FALSE, size = 1.5) +
        scale_x_continuous(trans = xTrans) +
        scale_y_continuous(trans = yTrans)

    # Label interesting groups
    p <- p + labs(color = paste(interestingGroups, collapse = ":\n"))

    # Color palette
    if (is(color, "ScaleDiscrete")) {
        p <- p + color
    }

    # Facets
    facets <- NULL
    if (.isAggregate(object)) {
        facets <- c(facets, "sampleNameAggregate")
    }
    if (is.character(facets)) {
        p <- p + facet_wrap(facets = facets, scales = "free")
    }

    p
}
