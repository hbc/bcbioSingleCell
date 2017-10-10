#' Plot Mitochondrial Transcript Abundance
#'
#' @rdname plotMitoRatio
#' @name plotMitoRatio
#' @family Quality Control Metrics
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @inherit plotGenesPerCell
NULL



# Constructors ====
.plotMitoRatioBoxplot <- function(
    object,
    interestingGroup = "sampleName",
    max = 1,
    filterCells = FALSE,
    aggregateReplicates = TRUE) {
    metrics <- metrics(
        object,
        filterCells = filterCells,
        aggregateReplicates = aggregateReplicates)
    p <- ggplot(
        metrics,
        mapping = aes_string(
            x = "sampleName",
            y = "mitoRatio",
            fill = interestingGroup)
    ) +
        scale_y_sqrt() +
        labs(x = "sample",
             y = "mito ratio") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))

    if (!isTRUE(aggregateReplicates) &
        "sampleNameAggregate" %in% colnames(metrics) &
        interestingGroup == "sampleName") {
        p <- p +
            geom_boxplot(
                color = lineColor,
                fill = "white")
    } else {
        p <- p +
            geom_boxplot(
                alpha = qcPlotAlpha,
                color = lineColor) +
            scale_fill_viridis(discrete = TRUE)
    }

    # Median labels
    if (length(unique(metrics[["sampleName"]])) <= qcLabelMaxNum) {
        formula <- formula(paste("mitoRatio", "sampleName", sep = " ~ "))
        meta <- sampleMetadata(
            object,
            aggregateReplicates = aggregateReplicates)
        medianMitoRatio <-
            aggregate(
                formula = formula,
                data = metrics,
                FUN = median) %>%
            left_join(meta, by = "sampleName") %>%
            mutate(mitoRatio = round(.data[["mitoRatio"]], digits = 3))
        p <- p +
            geom_label(
                data = medianMitoRatio,
                mapping = aes_string(label = "mitoRatio"),
                alpha = qcLabelAlpha,
                color = qcLabelColor,
                fill = qcLabelFill,
                fontface = qcLabelFontface,
                label.padding = qcLabelPadding,
                label.size = qcLabelSize,
                show.legend = FALSE)
    }

    # Cutoff lines
    if (max < 1) {
        p <- p +
            .qcCutoffLine(yintercept = max)
    }

    # Facets
    facets <- NULL
    if (isTRUE(metadata(object)[["multiplexedFASTQ"]])) {
        facets <- c(facets, "description")
    }
    if (!isTRUE(aggregateReplicates) &
        "sampleNameAggregate" %in% colnames(metrics)) {
        facets <- c(facets, "sampleNameAggregate")
        if (interestingGroup == "sampleName") {
            p <- p +
                theme(legend.position = "none")
        }
    }
    if (!is.null(facets)) {
        p <- p +
            facet_wrap(facets = facets,
                       scales = "free_x")
    }

    p
}



.plotMitoRatioRidgeline <- function(
    object,
    interestingGroup = "sampleName",
    max = 1,
    filterCells = FALSE,
    aggregateReplicates = TRUE) {
    metrics <- metrics(
        object,
        filterCells = filterCells,
        aggregateReplicates = aggregateReplicates)
    p <- ggplot(
        metrics,
        mapping = aes_string(
            x = "mitoRatio",
            y = "sampleName",
            fill = interestingGroup)
    ) +
        labs(x = "mito ratio",
             y = "sample") +
        geom_density_ridges(
            alpha = qcPlotAlpha,
            color = lineColor,
            panel_scaling = TRUE,
            scale = qcRidgeScale) +
        scale_fill_viridis(discrete = TRUE) +
        scale_x_sqrt() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))

    # Cutoff lines
    if (max < 1) {
        p <- p +
            .qcCutoffLine(xintercept = max)
    }

    # Facets
    facets <- NULL
    if (isTRUE(metadata(object)[["multiplexedFASTQ"]])) {
        facets <- c(facets, "description")
    }
    if (!isTRUE(aggregateReplicates) &
        "sampleNameAggregate" %in% colnames(metrics)) {
        facets <- c(facets, "sampleNameAggregate")
        # Turn off the legend
        p <- p +
            theme(legend.position = "none")
    }
    if (!is.null(facets)) {
        p <- p +
            facet_wrap(facets = facets)
    }

    p
}



.plotMitoRatioScatterplot <- function(
    object,
    interestingGroup = "sampleName",
    filterCells = FALSE,
    aggregateReplicates = TRUE) {
    metrics <- metrics(
        object,
        filterCells = filterCells,
        aggregateReplicates = aggregateReplicates)
    p <- ggplot(
        metrics,
        mapping = aes_string(
            x = "nCoding",
            y = "nMito",
            color = interestingGroup)
    ) +
        labs(x = "mito counts",
             y = "coding counts") +
        scale_x_sqrt() +
        scale_y_sqrt()

    if (!isTRUE(aggregateReplicates) &
        "sampleNameAggregate" %in% colnames(metrics) &
        interestingGroup == "sampleName") {
        p <- p +
            geom_point(
                color = "gray",
                size = 0.8) +
            geom_smooth(
                color = "black",
                method = "gam",
                se = FALSE,
                size = 1.5)
    } else {
        p <- p +
            geom_point(
                alpha = 0.25,
                size = 0.8) +
            geom_smooth(
                method = "gam",
                se = FALSE,
                size = 1.5) +
            scale_color_viridis(discrete = TRUE)
    }

    # Facets
    facets <- NULL
    if (isTRUE(metadata(object)[["multiplexedFASTQ"]])) {
        facets <- c(facets, "description")
    }
    if (!isTRUE(aggregateReplicates) &
        "sampleNameAggregate" %in% colnames(metrics)) {
        facets <- c(facets, "sampleNameAggregate")
        # Turn off the legend
        p <- p +
            theme(legend.position = "none")
    }
    if (!is.null(facets)) {
        p <- p +
            facet_wrap(facets = facets)
    }

    p
}



.plotMitoRatio <- function(
    object,
    interestingGroup,
    max,
    filterCells = FALSE,
    aggregateReplicates = TRUE) {
    if (missing(interestingGroup)) {
        interestingGroup <- interestingGroups(object)[[1]]
    }
    if (missing(max)) {
        max <- object %>%
            metadata() %>%
            .[["filterParams"]] %>%
            .[["maxMitoRatio"]]
        if (is.null(max)) {
            max <- 1
        }
    }
    suppressMessages(plot_grid(
        .plotMitoRatioRidgeline(
            object,
            interestingGroup = interestingGroup,
            max = max,
            filterCells = filterCells,
            aggregateReplicates = aggregateReplicates),
        .plotMitoRatioBoxplot(
            object,
            interestingGroup = interestingGroup,
            max = max,
            filterCells = filterCells,
            aggregateReplicates = aggregateReplicates),
        .plotMitoRatioScatterplot(
            object,
            filterCells = filterCells,
            aggregateReplicates = aggregateReplicates),
        labels = "auto",
        nrow = 3
    ))
}



# Methods ====
#' @rdname plotMitoRatio
#' @export
setMethod(
    "plotMitoRatio",
    signature("bcbioSingleCellANY"),
    .plotMitoRatio)
