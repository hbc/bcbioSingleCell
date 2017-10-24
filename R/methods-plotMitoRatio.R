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
#' @importFrom viridis scale_fill_viridis
.plotMitoRatioBoxplot <- function(
    object,
    interestingGroups = "sampleName",
    max = 1,
    filterCells = TRUE,
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
            fill = interestingGroups)
    ) +
        scale_y_sqrt() +
        labs(x = "sample",
             y = "mito ratio") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))

    if (!isTRUE(aggregateReplicates) &
        "sampleNameAggregate" %in% colnames(metrics) &
        interestingGroups == "sampleName") {
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
            left_join(meta, by = "sampleName")
        p <- p +
            geom_label(
                data = medianMitoRatio,
                mapping = aes_(label = ~round(mitoRatio, digits = 2)),
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
    if (isTRUE(metadata(object)[["multiplexedFASTQ"]]) &
        length(unique(metrics[["description"]])) > 1) {
        facets <- c(facets, "description")
    }
    if (!isTRUE(aggregateReplicates) &
        "sampleNameAggregate" %in% colnames(metrics)) {
        facets <- c(facets, "sampleNameAggregate")
        if (interestingGroups == "sampleName") {
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



#' @importFrom ggridges geom_density_ridges
.plotMitoRatioRidgeline <- function(
    object,
    interestingGroups = "sampleName",
    max = 1,
    filterCells = TRUE,
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
            fill = interestingGroups)
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
    if (isTRUE(metadata(object)[["multiplexedFASTQ"]]) &
        length(unique(metrics[["description"]])) > 1) {
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



#' @importFrom viridis scale_color_viridis
.plotMitoRatioScatterplot <- function(
    object,
    interestingGroups = "sampleName",
    filterCells = TRUE,
    aggregateReplicates = TRUE) {
    metrics <- metrics(
        object,
        filterCells = filterCells,
        aggregateReplicates = aggregateReplicates) %>%
        rownames_to_column("cellID")
        # Drop cells with zeroes for both `nCoding` and `nMito`
    dropCells <- dplyr::filter(
        metrics,
        .data[["nCoding"]] == 0 &
        .data[["nMito"]] == 0) %>%
        pull("cellID")
    metrics <- dplyr::filter(metrics, !.data[["cellID"]] %in% dropCells)
    if (nrow(metrics) == 0) {
        stop(paste(
            "No cells contain coding and mito counts.",
            "Check that your organism is set correctly",
            "and rerun 'calculateMetrics()'."
            ), call. = FALSE)
    }
    p <- ggplot(
        metrics,
        mapping = aes_string(
            x = "nCoding",
            y = "nMito",
            color = interestingGroups)
    ) +
        labs(x = "mito counts",
             y = "coding counts") +
        scale_x_sqrt() +
        scale_y_sqrt()

    if (!isTRUE(aggregateReplicates) &
        "sampleNameAggregate" %in% colnames(metrics) &
        interestingGroups == "sampleName") {
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
    if (isTRUE(metadata(object)[["multiplexedFASTQ"]]) &
        length(unique(metrics[["description"]])) > 1) {
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



#' @importFrom cowplot plot_grid
.plotMitoRatio <- function(
    object,
    interestingGroups,
    max,
    filterCells = TRUE,
    aggregateReplicates = TRUE) {
    if (missing(interestingGroups)) {
        interestingGroups <-
            metadata(object)[["interestingGroups"]][[1]]
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
            interestingGroups = interestingGroups,
            max = max,
            filterCells = filterCells,
            aggregateReplicates = aggregateReplicates),
        .plotMitoRatioBoxplot(
            object,
            interestingGroups = interestingGroups,
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
    signature("bcbioSingleCell"),
    .plotMitoRatio)
