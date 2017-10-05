#' Plot UMIs per Cell
#'
#' Plot the universal molecular identifiers (UMIs) per cell.
#'
#' @rdname plotUMIsPerCell
#' @name plotUMIsPerCell
#' @family Quality Control Metrics
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @inherit plotGenesPerCell
NULL



# Constructors ====
.plotUMIsPerCellBoxplot <- function(
    object,
    interestingGroup = "sampleName",
    min = NULL,
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
            y = "nUMI",
            fill = interestingGroup)
    ) +
        labs(x = "sample",
             y = "umis per cell") +
        geom_boxplot(color = lineColor) +
        scale_y_log10() +
        scale_fill_viridis(discrete = TRUE) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))

    # Median labels
    if (length(unique(metrics[["sampleName"]])) <= qcLabelMaxNum) {
        formula <- formula(paste("nUMI", "sampleName", sep = " ~ "))
        meta <- sampleMetadata(
            object,
            aggregateReplicates = aggregateReplicates)
        medianUMIs <- aggregate(
            formula = formula,
            data = metrics,
            FUN = median) %>%
            left_join(meta, by = "sampleName") %>%
            mutate(nUMI = round(.data[["nUMI"]]))
        p <- p +
            geom_label(
                data = medianUMIs,
                mapping = aes_string(label = "nUMI"),
                alpha = qcLabelAlpha,
                color = qcLabelColor,
                fill = qcLabelFill,
                fontface = qcLabelFontface,
                label.padding = qcLabelPadding,
                label.size = qcLabelSize,
                show.legend = FALSE)
    }

    # Cutoff lines
    if (!is.null(min)) {
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



.plotUMIsPerCellHistogram <- function(
    object,
    min = NULL,
    filterCells = FALSE,
    aggregateReplicates = FALSE) {
    metrics <- metrics(
        object,
        filterCells = filterCells,
        aggregateReplicates = aggregateReplicates)
    p <- ggplot(
        metrics,
        mapping = aes_string(
            x = "nUMI",
            fill = "sampleName")
    ) +
        labs(x = "umis per cell",
             fill = "sample") +
        geom_histogram(bins = bins) +
        scale_x_log10() +
        scale_y_sqrt() +
        scale_fill_viridis(discrete = TRUE) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))

    # Cutoff lines
    if (!is.null(min)) {
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



.plotUMIsPerCell <- function(
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
            .[["minUMIs"]]
    }
    plot_grid(
        .plotUMIsPerCellHistogram(
            object,
            min = min,
            filterCells = filterCells,
            aggregateReplicates = aggregateReplicates),
        .plotUMIsPerCellBoxplot(
            object,
            interestingGroup = interestingGroup,
            min = min,
            filterCells = filterCells,
            aggregateReplicates = aggregateReplicates),
        labels = "auto",
        nrow = 2
    )
}



# Methods ====
#' @rdname plotUMIsPerCell
#' @export
setMethod("plotUMIsPerCell", "bcbioSingleCell", .plotUMIsPerCell)
