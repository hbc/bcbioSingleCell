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
    p <- ggplot(cellCounts,
                aes_(x = ~sampleName,
                     y = ~cells,
                     fill = as.name(interestingGroup))) +
        labs(x = "sample",
             y = "cell count") +
        geom_bar(stat = "identity") +
        geom_text(
            aes_(label = ~cells),
            fontface = "bold",
            vjust = -0.5) +
        scale_fill_viridis(discrete = TRUE) +
        theme(axis.text.x = element_text(angle = 15, hjust = 1))
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
