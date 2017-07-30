#' Plot Cell Counts
#'
#' @rdname plotCellCounts
#' @name plotCellCounts
#' @author Michael Steinbaugh, Rory Kirchner
#' @inherit plotGenesPerCell
NULL



# Constructors ====
.plotCellCounts <- function(object) {
    metrics <- metrics(object)
    cellCounts <- metrics %>%
        group_by(!!sym("sampleID")) %>%
        summarise(cells = n()) %>%
        left_join(sampleMetadata(object), by = "sampleID")
    interestingGroup <- interestingGroups(object)[[1L]]
    p <- ggplot(cellCounts,
                aes_(x = ~sampleName,
                     y = ~cells,
                     fill = as.name(interestingGroup))) +
        labs(x = "sample",
             y = "cell count") +
        geom_bar(stat = "identity") +
        geom_text(vjust = -0.5, aes_(label = ~cells)) +
        theme(axis.text.x = element_text(angle = 90L, hjust = 1L))
    if (isTRUE(metadata(object)[["multiplexedFASTQ"]])) {
        p <- p + facet_wrap(~fileName)
    }
    p
}



# Methods ====
#' @rdname plotCellCounts
#' @export
setMethod("plotCellCounts", "bcbioSCDataSet", .plotCellCounts)



#' @rdname plotCellCounts
#' @export
setMethod("plotCellCounts", "bcbioSCSubset", .plotCellCounts)
