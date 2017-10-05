#' Plot UMI and Gene Correlation
#'
#' @rdname plotUMIsVsGenes
#' @name plotUMIsVsGenes
#' @family Quality Control Metrics
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @inherit plotGenesPerCell
NULL



# Constructors ====
.plotUMIsVsGenes <- function(
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
            x = "nUMI",
            y = "nGene",
            color = interestingGroup)
    ) +
        labs(x = "umis per cell",
             y = "genes per cell") +
        geom_point(alpha = 0.25, size = 0.8) +
        geom_smooth(method = "gam", se = FALSE, size = 2) +
        scale_x_log10() +
        scale_y_log10() +
        scale_color_viridis(discrete = TRUE) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))

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



# Methods ====
#' @rdname plotUMIsVsGenes
#' @export
setMethod("plotUMIsVsGenes", "bcbioSingleCellANY", .plotUMIsVsGenes)
