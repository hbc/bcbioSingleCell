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
    aggregateReplicates = FALSE) {
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
        geom_boxplot(color = lineColor) +
        scale_y_sqrt() +
        scale_fill_viridis(discrete = TRUE) +
        labs(x = "sample",
             y = "mito ratio") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
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
    if (isTRUE(metadata(object)[["multiplexedFASTQ"]])) {
        facets <- c(facets, "fileName")
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



.plotMitoRatioHistogram <- function(
    object,
    interestingGroup = "sampleName",
    max = 1,
    filterCells = FALSE,
    aggregateReplicates = FALSE) {
    metrics <- metrics(
        object,
        filterCells = filterCells,
        aggregateReplicates = aggregateReplicates)
    p <- ggplot(
        metrics,
        mapping = aes_string(
            x = "mitoRatio",
            fill = interestingGroup)
    ) +
        geom_histogram(bins = bins) +
        scale_x_sqrt() +
        scale_y_sqrt() +
        scale_fill_viridis(discrete = TRUE) +
        labs(x = "mito ratio",
             y = "sample")

    # Cutoff lines
    if (max < 1) {
        p <- p +
            .qcCutoffLine(xintercept = max)
    }

    # Facets
    facets <- NULL
    if (isTRUE(metadata(object)[["multiplexedFASTQ"]])) {
        facets <- c(facets, "fileName")
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
    aggregateReplicates = FALSE) {
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
        geom_point(alpha = 0.25, size = 0.8) +
        geom_smooth(method = "gam", se = FALSE, size = 2) +
        scale_x_sqrt() +
        scale_y_sqrt() +
        scale_color_viridis(discrete = TRUE)

    # Facets
    facets <- NULL
    if (isTRUE(metadata(object)[["multiplexedFASTQ"]])) {
        facets <- c(facets, "fileName")
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
    aggregateReplicates = FALSE) {
    if (missing(interestingGroup)) {
        interestingGroup <- interestingGroups(object)[[1]]
    }
    if (missing(max)) {
        max <- object %>%
            metadata() %>%
            .[["filterParams"]] %>%
            .[["maxMitoRatio"]]
    }
    plot_grid(
        .plotMitoRatioScatterplot(
            object,
            filterCells = filterCells,
            aggregateReplicates = aggregateReplicates),
        .plotMitoRatioHistogram(
            object,
            max = max,
            filterCells = filterCells,
            aggregateReplicates = aggregateReplicates),
        .plotMitoRatioBoxplot(
            object,
            interestingGroup = interestingGroup,
            max = max,
            filterCells = filterCells,
            aggregateReplicates = aggregateReplicates),
        labels = "auto",
        nrow = 3)
}



# Methods ====
#' @rdname plotMitoRatio
#' @export
setMethod("plotMitoRatio", "bcbioSingleCellANY", .plotMitoRatio)
