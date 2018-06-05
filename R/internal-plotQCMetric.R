# Plot a single quality control metric
.plotQCMetric <- function(
    object,
    metricCol,
    geom = c("ecdf", "ridgeline", "violin", "histogram", "boxplot"),
    interestingGroups,
    min = 0L,
    max = Inf,
    trans = "identity",
    ratio = FALSE,
    color = scale_color_hue(),
    fill = scale_fill_hue(),
    title = NULL
) {
    assert_is_a_string(metricCol)
    geom <- match.arg(geom)
    if (missing(interestingGroups)) {
        interestingGroups <- bcbioBase::interestingGroups(object)
    } else {
        interestingGroups(object) <- interestingGroups
    }
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
    trendline = FALSE,
    color = scale_color_hue(),
    title = NULL
) {
    assert_is_a_string(xCol)
    assert_is_a_string(yCol)
    assert_is_a_string(xTrans)
    assert_is_a_string(yTrans)
    if (missing(interestingGroups)) {
        interestingGroups <- bcbioBase::interestingGroups(object)
    } else {
        interestingGroups(object) <- interestingGroups
    }
    assertIsColorScaleDiscreteOrNULL(color)
    assertIsAStringOrNULL(title)

    data <- metrics(
        object = object,
        interestingGroups = interestingGroups,
        prefilter = FALSE
    )
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
        mapping = aes_string(
            x = xCol,
            y = yCol,
            color = "interestingGroups"
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
        p <- p + facet_wrap(facets = facets, scales = "free")
    }

    p
}
