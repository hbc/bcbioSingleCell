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
    aggregateReplicates = TRUE) {
    metrics <- metrics(
        object,
        filterCells = filterCells,
        aggregateReplicates = aggregateReplicates)
    p <- ggplot(
        metrics,
        mapping = aes_string(
            x = "nUMI",
            y = "nGene",
            color = interestingGroup,
            fill = interestingGroup)
    ) +
        labs(x = "umis per cell",
             y = "genes per cell") +
        scale_x_log10() +
        scale_y_log10() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))

    if (!isTRUE(aggregateReplicates) &
        "sampleNameAggregate" %in% colnames(metrics) &
        interestingGroup == "sampleName") {
        p <- p +
            geom_point(
                color = "gray",
                size = 0.8) +
            geom_smooth(
                color = "black",
                method = "gam",
                se = FALSE,
                size = 1.5)
    } else {
        p <- p +
            geom_point(
                alpha = 0.25,
                size = 0.8) +
            geom_smooth(
                method = "gam",
                se = FALSE,
                size = 1.5) +
            scale_color_viridis(discrete = TRUE)
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



# Methods ====
#' @rdname plotUMIsVsGenes
#' @export
setMethod("plotUMIsVsGenes", "bcbioSingleCellANY", .plotUMIsVsGenes)
