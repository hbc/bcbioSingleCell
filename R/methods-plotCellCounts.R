#' Plot Cell Counts
#'
#' @rdname plotCellCounts
#' @name plotCellCounts
#' @family Quality Control Metrics
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @inherit plotGenesPerCell
NULL



# Constructors ====
.plotCellCounts <- function(
    object,
    interestingGroup,
    filterCells = FALSE,
    aggregateReplicates = TRUE) {
    if (missing(interestingGroup)) {
        interestingGroup <- .interestingGroup(object)
    }
    metrics <- metrics(
        object,
        filterCells = filterCells,
        aggregateReplicates = aggregateReplicates)
    meta <- sampleMetadata(
        object,
        aggregateReplicates = aggregateReplicates)
    cellCounts <- metrics %>%
        group_by(!!sym("sampleName")) %>%
        summarize(cells = n()) %>%
        left_join(meta, by = "sampleName")
    p <- ggplot(
        cellCounts,
        mapping = aes_string(
            x = "sampleName",
            y = "cells",
            fill = interestingGroup)
    ) +
        labs(x = "sample",
             y = "cell count") +
        geom_bar(stat = "identity") +
        scale_fill_viridis(discrete = TRUE) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))

    # Labels
    if (nrow(cellCounts) <= qcLabelMaxNum) {
        p <- p +
            geom_text(
                mapping = aes_string(label = "cells"),
                fontface = "bold",
                vjust = -0.5)
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
    }
    if (!is.null(facets)) {
        p <- p +
            facet_wrap(facets = facets,
                       scales = "free_x")
    }

    p
}



# Methods ====
#' @rdname plotCellCounts
#' @export
setMethod("plotCellCounts", "bcbioSingleCellANY", .plotCellCounts)
