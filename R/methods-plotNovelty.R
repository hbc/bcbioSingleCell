#' Plot Novelty Score
#'
#' @rdname plotNovelty
#' @name plotNovelty
#' @family Quality Control Metrics
#' @author Michael Steinbaugh
#'
#' @inherit plotGenesPerCell
#'
#' @details "Novelty" refers to log10 genes detected per count.
NULL



# Constructors ====
.plotNoveltyBoxplot <- function(
    object,
    interestingGroup = "sampleName",
    min = 0,
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
            y = "log10GenesPerUMI",
            fill = interestingGroup)
    ) +
        labs(x = "sample",
             y = "log10 genes per UMI") +
        scale_y_sqrt() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))

    if (!isTRUE(aggregateReplicates) &
        "sampleNameAggregate" %in% colnames(metrics) &
        interestingGroup == "sampleName") {
        p <- p +
            geom_boxplot(color = lineColor, fill = "white")
    } else {
        p <- p +
            geom_boxplot(color = lineColor) +
            scale_fill_viridis(discrete = TRUE)
    }

    # Median labels
    if (length(unique(metrics[["sampleName"]])) <= qcLabelMaxNum) {
        formula <- formula(paste("log10GenesPerUMI", "sampleName", sep = " ~ "))
        meta <- sampleMetadata(
            object,
            aggregateReplicates = aggregateReplicates)
        medianNovelty <-
            aggregate(
                formula = formula,
                data = metrics,
                FUN = median) %>%
            left_join(meta, by = "sampleName")
        p <- p +
            geom_label(
                data = medianNovelty,
                aes_(label = ~round(log10GenesPerUMI, digits = 2)),
                alpha = qcLabelAlpha,
                color = qcLabelColor,
                fill = qcLabelFill,
                fontface = qcLabelFontface,
                label.padding = qcLabelPadding,
                label.size = qcLabelSize,
                show.legend = FALSE)
    }

    # Cutoff lines
    if (min > 0) {
        p <- p +
            .qcCutoffLine(yintercept = min)
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



.plotNoveltyHistogram <- function(
    object,
    min = 0,
    filterCells = FALSE,
    aggregateReplicates = FALSE) {
    metrics <- metrics(
        object,
        filterCells = filterCells,
        aggregateReplicates = aggregateReplicates)
    p <- ggplot(
        metrics,
        mapping = aes_string(
            x = "log10GenesPerUMI",
            fill = "sampleName")
    ) +
        labs(x = "log10 genes per UMI",
             fill = "sample") +
        geom_histogram(bins = bins) +
        scale_x_sqrt() +
        scale_y_sqrt() +
        scale_fill_viridis(discrete = TRUE)

    # Cutoff lines
    if (min > 0) {
        p <- p +
            .qcCutoffLine(xintercept = min)
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



.plotNovelty <- function(
    object,
    interestingGroup,
    min,
    filterCells = FALSE,
    aggregateReplicates = FALSE) {
    if (missing(interestingGroup)) {
        interestingGroup <- interestingGroups(object)[[1]]
    }
    if (missing(min)) {
        min <- object %>%
            metadata() %>%
            .[["filterParams"]] %>%
            .[["minNovelty"]]
    }
    plot_grid(
        .plotNoveltyHistogram(
            object,
            min = min,
            filterCells = filterCells,
            aggregateReplicates = aggregateReplicates),
        .plotNoveltyBoxplot(
            object,
            interestingGroup = interestingGroup,
            min = min,
            filterCells = filterCells,
            aggregateReplicates = aggregateReplicates),
        labels = "auto",
        nrow = 2)
}



# Methods ====
#' @rdname plotNovelty
#' @export
setMethod("plotNovelty", "bcbioSingleCellANY", .plotNovelty)
