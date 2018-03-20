# FIXME PBMC working example is broken

#' Plot Cell Counts
#'
#' @name plotCellCounts
#' @family Quality Control Functions
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @inherit plotGenesPerCell
#'
#' @param metadata Sample metadata `data.frame`.
#'
#' @examples
#' # bcbioSingleCell ====
#' plotCellCounts(bcb_small)
#'
#' # seurat ====
#' plotCellCounts(seurat_small)
#' plotCellCounts(pbmc_small)
NULL



# Constructors =================================================================
#' @importFrom bcbioBase interestingGroups
#' @importFrom dplyr group_by left_join n summarize
#' @importFrom ggplot2 aes_string element_text facet_wrap geom_bar geom_label
#'   ggplot labs theme
#' @importFrom rlang !! sym
#' @importFrom viridis scale_fill_viridis
.plotCellCounts <- function(
    object,
    interestingGroups,
    fill = scale_fill_viridis(discrete = TRUE)
) {
    if (missing(interestingGroups)) {
        interestingGroups <- bcbioBase::interestingGroups(object)
    }

    # sampleID level issue here
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
        theme(axis.text.x = element_text(angle = 90L, hjust = 1L))

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
    if (isTRUE(.checkAggregate(data))) {
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
