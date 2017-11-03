#' Plot Cell Counts
#'
#' @rdname plotCellCounts
#' @name plotCellCounts
#' @family Quality Control Metrics
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @inherit plotGenesPerCell
#'
#' @examples
#' # bcbioSingleCell
#' \dontrun{
#' plotCellCounts(bcb)
#' }
#'
#' # seurat
#' \dontrun{
#' plotCellCounts(seurat)
#' }
#'
#' # metrics data.frame
#' \dontrun{
#' metrics <- metrics(bcb)
#' metadata <- sampleMetadata(bcb)
#' plotCellCounts(metrics, metadata)
#' }
NULL



# Constructors ====
#' @importFrom dplyr group_by left_join n summarize
#' @importFrom ggplot2 aes_string element_text facet_wrap geom_bar geom_text
#'   ggplot labs theme
#' @importFrom rlang !! sym
#' @importFrom viridis scale_fill_viridis
.plotCellCounts <- function(
    object,
    metadata,
    interestingGroups,
    multiplexed = FALSE,
    fill = scale_fill_viridis(discrete = TRUE)) {
    data <- object %>%
        group_by(!!sym("sampleName")) %>%
        summarize(nCells = n()) %>%
        left_join(metadata, by = "sampleName")
    p <- ggplot(
        data,
        mapping = aes_string(
            x = "sampleName",
            y = "nCells",
            fill = "interestingGroups")
    ) +
        geom_bar(stat = "identity") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))

    # Label interesting groups
    if (!missing(interestingGroups)) {
        p <- p + labs(fill = paste(interestingGroups, collapse = ":\n"))
    } else {
        p <- p + labs(color = NULL, fill = NULL)
    }

    # Color palette
    if (!is.null(fill)) {
        p <- p + fill
    }

    # Labels
    if (nrow(data) <= qcLabelMaxNum) {
        p <- p +
            geom_text(
                mapping = aes_string(label = "nCells"),
                fontface = "bold",
                vjust = -0.5)
    }

    # Facets
    facets <- NULL
    if (isTRUE(multiplexed) &
        length(unique(object[["description"]])) > 1) {
        facets <- c(facets, "description")
    }
    if (!isTRUE(aggregateReplicates) &
        "sampleNameAggregate" %in% colnames(object)) {
        facets <- c(facets, "sampleNameAggregate")
    }
    if (!is.null(facets)) {
        p <- p + facet_wrap(facets = facets, scales = "free_x")
    }

    p
}



# Methods ====
#' @rdname plotCellCounts
#' @export
setMethod(
    "plotCellCounts",
    signature("bcbioSingleCell"),
    function(
        object,
        interestingGroups,
        filterCells = TRUE,
        aggregateReplicates = TRUE,
        fill = scale_fill_viridis(discrete = TRUE)) {
        if (missing(interestingGroups)) {
            interestingGroups <-
                metadata(object)[["interestingGroups"]]
        }
        metrics <- metrics(
            object,
            interestingGroups = interestingGroups,
            filterCells = filterCells,
            aggregateReplicates = aggregateReplicates)
        metadata <- sampleMetadata(
            object,
            aggregateReplicates = aggregateReplicates)
        multiplexed <- metadata(object)[["multiplexedFASTQ"]]
        .plotCellCounts(
            object = metrics,
            metadata = metadata,
            interestingGroups = interestingGroups,
            multiplexed = multiplexed,
            fill = fill)
    })



#' @rdname plotCellCounts
#' @export
setMethod(
    "plotCellCounts",
    signature("data.frame"),
    .plotCellCounts)



#' @rdname plotCellCounts
#' @export
setMethod(
    "plotCellCounts",
    signature("seurat"),
    function(
        object,
        interestingGroups,
        fill = scale_fill_viridis(discrete = TRUE)) {
        if (missing(interestingGroups)) {
            interestingGroups <-
                slot(object, "misc") %>%
                .[["bcbio"]] %>%
                .[["interestingGroups"]]
        }
        metrics <- metrics(object, interestingGroups = interestingGroups)
        metadata <- sampleMetadata(object)
        .plotCellCounts(
            object = metrics,
            metadata = metadata,
            interestingGroups = interestingGroups,
            multiplexed = FALSE,
            fill = fill)
    })
