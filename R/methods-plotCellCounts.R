#' Plot Cell Counts
#'
#' @name plotCellCounts
#' @family Quality Control Metrics
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @inheritParams general
#'
#' @return `ggplot`.
#'
#' @examples
#' load(system.file("extdata/bcb.rda", package = "bcbioSingleCell"))
#' load(system.file("extdata/seurat.rda", package = "bcbioSingleCell"))
#'
#' # bcbioSingleCell ====
#' plotCellCounts(bcb)
#'
#' # seurat ====
#' plotCellCounts(seurat)
NULL



# Constructors =================================================================
.plotCellCounts <- function(
    object,
    interestingGroups,
    fill = scale_fill_viridis(discrete = TRUE)
) {
    if (missing(interestingGroups)) {
        interestingGroups <- bcbioBase::interestingGroups(object)
    }
    assertIsFillScaleDiscreteOrNULL(fill)

    metrics <- metrics(object, interestingGroups)
    metadata <- sampleData(object, interestingGroups)

    data <- metrics %>%
        group_by(!!sym("sampleName")) %>%
        summarize(nCells = n()) %>%
        left_join(metadata, by = "sampleName")

    p <- ggplot(
        data = data,
        mapping = aes_string(
            x = "sampleName",
            y = "nCells",
            fill = "interestingGroups"
        )
    ) +
        geom_bar(stat = "identity") +
        labs(fill = paste(interestingGroups, collapse = ":\n")) +
        theme(axis.text.x = element_text(angle = 90L, hjust = 1L))

    # Color palette
    if (!is.null(fill)) {
        p <- p + fill
    }

    # Labels
    if (nrow(data) <= qcLabelMaxNum) {
        p <- p + geom_label(
            data = data,
            mapping = aes_string(label = "nCells"),
            alpha = qcLabelAlpha,
            color = qcLabelColor,
            fill = qcLabelFill,
            fontface = qcLabelFontface,
            label.padding = qcLabelPadding,
            label.size = qcLabelSize,
            show.legend = FALSE,
            # Align the label just under the top of the bar
            vjust = 1.25
        )
    }

    # Facets
    facets <- NULL
    if (isTRUE(.checkAggregate(object))) {
        facets <- c(facets, "sampleNameAggregate")
    }
    if (is.character(facets)) {
        p <- p + facet_wrap(facets = facets, scales = "free_x")
    }

    p
}



# Methods ======================================================================
#' @rdname plotCellCounts
#' @export
setMethod(
    "plotCellCounts",
    signature("bcbioSingleCell"),
    .plotCellCounts
)



#' @rdname plotCellCounts
#' @export
setMethod(
    "plotCellCounts",
    signature("seurat"),
    .plotCellCounts
)
