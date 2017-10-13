#' Plot Genes per Cell
#'
#' @rdname plotGenesPerCell
#' @name plotGenesPerCell
#' @family Quality Control Metrics
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @inheritParams AllGenerics
#' @inheritParams metrics
#' @param interestingGroups Interesting group, to use for colors.
#' @param min Recommended minimum value cutoff.
#' @param max Recommended maximum value cutoff.
#'
#' @return [ggplot] grid.
NULL



# Constructors ====
.plotGenesPerCellBoxplot <- function(
    object,
    interestingGroups = "sampleName",
    min = 0,
    max = Inf,
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
            y = "nGene",
            fill = interestingGroups)
    ) +
        labs(x = "sample",
             y = "genes per cell") +
        scale_y_sqrt() +
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
        formula <- formula(paste("nGene", "sampleName", sep = " ~ "))
        meta <- sampleMetadata(
            object,
            aggregateReplicates = aggregateReplicates)
        medianGenes <-
            aggregate(
                formula = formula,
                data = metrics,
                FUN = median) %>%
            left_join(meta, by = "sampleName")
        p <- p +
            geom_label(
                data = medianGenes,
                mapping = aes_(label = ~round(nGene)),
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
    if (max < Inf) {
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



.plotGenesPerCellRidgeline <- function(
    object,
    interestingGroups = "sampleName",
    min = 0,
    max = Inf,
    filterCells = TRUE,
    aggregateReplicates = TRUE) {
    metrics <- metrics(
        object,
        filterCells = filterCells,
        aggregateReplicates = aggregateReplicates)
    p <- ggplot(
        metrics,
        mapping = aes_string(
            x = "nGene",
            y = "sampleName",
            fill = interestingGroups)
    ) +
        labs(x = "genes per cell",
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
    if (max < Inf) {
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
            facet_wrap(facets = facets,
                       scales = "free_y")
    }

    p
}



.plotGenesPerCell <- function(
    object,
    interestingGroups,
    min,
    max,
    filterCells = TRUE,
    aggregateReplicates = TRUE) {
    if (missing(interestingGroups)) {
        interestingGroups <-
            metadata(object)[["interestingGroups"]][[1]]
    }
    if (missing(min)) {
        min <- object %>%
            metadata() %>%
            .[["filterParams"]] %>%
            .[["minGenes"]]
        if (is.null(min)) {
            min <- 0
        }
    }
    if (missing(max)) {
        max <- object %>%
            metadata() %>%
            .[["filterParams"]] %>%
            .[["maxGenes"]]
        if (is.null(max)) {
            max <- Inf
        }
    }
    suppressMessages(plot_grid(
        .plotGenesPerCellRidgeline(
            object,
            interestingGroups = interestingGroups,
            min = min,
            max = max,
            filterCells = filterCells,
            aggregateReplicates = aggregateReplicates),
        .plotGenesPerCellBoxplot(
            object,
            interestingGroups = interestingGroups,
            min = min,
            max = max,
            filterCells = filterCells,
            aggregateReplicates = aggregateReplicates),
        labels = "auto",
        nrow = 2
    ))
}



# Methods ====
#' @rdname plotGenesPerCell
#' @export
setMethod(
    "plotGenesPerCell",
    signature("bcbioSingleCellANY"),
    .plotGenesPerCell)
