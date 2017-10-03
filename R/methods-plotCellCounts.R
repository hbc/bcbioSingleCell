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
    filterCells = FALSE) {
    if (missing(interestingGroup)) {
        interestingGroup <- interestingGroups(object)[[1]]
    }
    cellCounts <- metrics(object, filterCells = filterCells) %>%
        group_by(!!sym("sampleID")) %>%
        summarize(cells = n()) %>%
        left_join(sampleMetadata(object), by = "sampleID")
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
    # Show cell counts for up to 20 samples (too busy otherwise)
    if (nrow(cellCounts) <= 20) {
        p <- p +
            geom_text(
                mapping = aes_string(label = "cells"),
                fontface = "bold",
                vjust = -0.5)
    }
    if (isTRUE(metadata(object)[["multiplexedFASTQ"]])) {
        p <- p +
            facet_wrap(~fileName)
    }
    p
}



# Methods ====
#' @rdname plotCellCounts
#' @export
setMethod("plotCellCounts", "bcbioSingleCellANY", .plotCellCounts)
