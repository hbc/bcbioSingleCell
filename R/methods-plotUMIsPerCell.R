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
    min = 0,
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
            y = "nUMI",
            fill = interestingGroup)
    ) +
        labs(x = "sample",
             y = "umis per cell") +
        scale_y_sqrt() +
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
    if (min > 0) {
        p <- p +
            .qcCutoffLine(yintercept = min)
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



.plotUMIsPerCellRidgeline <- function(
    object,
    interestingGroup = "sampleName",
    min = 0,
    filterCells = FALSE,
    aggregateReplicates = TRUE) {
    metrics <- metrics(
        object,
        filterCells = filterCells,
        aggregateReplicates = aggregateReplicates)
    p <- ggplot(
        metrics,
        mapping = aes_string(
            x = "nUMI",
            y = "sampleName",
            fill = interestingGroup)
    ) +
        labs(x = "umis per cell",
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
    if (min > 0) {
        p <- p +
            .qcCutoffLine(xintercept = min)
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



.plotUMIsPerCell <- function(
    object,
    interestingGroup,
    min,
    filterCells = FALSE,
    aggregateReplicates = TRUE) {
    if (missing(interestingGroup)) {
        interestingGroup <- interestingGroups(object)[[1]]
    }
    if (missing(min)) {
        min <- object %>%
            metadata() %>%
            .[["filterParams"]] %>%
            .[["minUMIs"]]
        if (is.null(min)) {
            min <- 0
        }
    }
    suppressMessages(plot_grid(
        .plotUMIsPerCellRidgeline(
            object,
            interestingGroup = interestingGroup,
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
    ))
}



# Methods ====
#' @rdname plotUMIsPerCell
#' @export
setMethod(
    "plotUMIsPerCell",
    signature("bcbioSingleCell"),
    .plotUMIsPerCell)
