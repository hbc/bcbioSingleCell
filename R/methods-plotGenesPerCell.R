#' Plot Genes per Cell
#'
#' @rdname plotGenesPerCell
#' @name plotGenesPerCell
#' @family Quality Control Metrics
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @inheritParams AllGenerics
#' @inheritParams metrics
#' @param interestingGroup Interesting group, to use for colors.
#' @param min Recommended minimum value cutoff.
#' @param max Recommended maximum value cutoff.
#'
#' @return [ggplot] grid.
NULL



# Constructors ====
.plotGenesPerCellBoxplot <- function(
    object,
    interestingGroup = "sampleName",
    min = NULL,
    max = NULL,
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
            y = "nGene",
            fill = interestingGroup)
    ) +
        labs(x = "sample",
             y = "genes per cell") +
        geom_boxplot(color = lineColor) +
        scale_y_sqrt() +
        scale_fill_viridis(discrete = TRUE) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))

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
            left_join(meta, by = "sampleName") %>%
            mutate(nGene = round(.data[["nGene"]]))
        p <- p +
            geom_label(
                data = medianGenes,
                mapping = aes_string(label = "nGene"),
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
    if (!is.null(max)) {
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



.plotGenesPerCellHistogram <- function(
    object,
    min = NULL,
    max = NULL,
    filterCells = FALSE,
    aggregateReplicates = FALSE) {
    metrics <- metrics(
        object,
        filterCells = filterCells,
        aggregateReplicates = aggregateReplicates)
    p <- ggplot(
        metrics,
        mapping = aes_string(
            x = "nGene",
            fill = "sampleName")
    ) +
        labs(x = "genes per cell",
             fill = "sample") +
        geom_histogram(bins = bins) +
        scale_x_sqrt() +
        scale_y_sqrt() +
        scale_fill_viridis(discrete = TRUE) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))

    # Cutoff lines
    if (!is.null(min)) {
        p <- p +
            .qcCutoffLine(xintercept = min)
    }
    if (!is.null(max)) {
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



.plotGenesPerCell <- function(
    object,
    interestingGroup,
    min,
    max,
    filterCells = FALSE,
    aggregateReplicates = FALSE) {
    if (missing(interestingGroup)) {
        interestingGroup <- interestingGroups(object)[[1]]
    }
    if (missing(min)) {
        min <- object %>%
            metadata() %>%
            .[["filterParams"]] %>%
            .[["minGenes"]]
    }
    if (missing(max)) {
        max <- object %>%
            metadata() %>%
            .[["filterParams"]] %>%
            .[["maxGenes"]]
    }
    plot_grid(
        .plotGenesPerCellHistogram(
            object,
            min = min,
            max = max,
            filterCells = filterCells,
            aggregateReplicates = aggregateReplicates),
        .plotGenesPerCellBoxplot(
            object,
            interestingGroup = interestingGroup,
            min = min,
            max = max,
            filterCells = filterCells,
            aggregateReplicates = aggregateReplicates),
        labels = "auto",
        nrow = 2)
}



# Methods ====
#' @rdname plotGenesPerCell
#' @export
setMethod("plotGenesPerCell", "bcbioSingleCellANY", .plotGenesPerCell)
