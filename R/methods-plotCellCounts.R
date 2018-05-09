#' Plot Cell Counts
#'
#' @name plotCellCounts
#' @family Quality Control Functions
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @inheritParams general
#'
#' @return `ggplot`.
#'
#' @examples
#' # bcbioSingleCell ====
#' plotCellCounts(indrops_small)
#'
#' # SingleCellExperiment ====
#' plotCellCounts(cellranger_small)
#'
#' # seurat ====
#' plotCellCounts(Seurat::pbmc_small)
NULL



# Methods ======================================================================
#' @rdname plotCellCounts
#' @export
setMethod(
    "plotCellCounts",
    signature("SingleCellExperiment"),
    function(
        object,
        interestingGroups,
        fill = scale_fill_hue(),
        title = "cell counts"
    ) {
        if (missing(interestingGroups)) {
            interestingGroups <- bcbioBase::interestingGroups(object)
        }
        assertIsFillScaleDiscreteOrNULL(fill)
        assertIsAStringOrNULL(title)

        metrics <- metrics(object, interestingGroups)
        sampleData <- sampleData(
            object = object,
            clean = FALSE,
            interestingGroups = interestingGroups,
            return = "data.frame"
        )
        sampleData[["sampleID"]] <- rownames(sampleData)
        data <- metrics %>%
            group_by(!!sym("sampleID")) %>%
            summarize(nCells = n()) %>%
            merge(sampleData, all.x = TRUE)

        p <- ggplot(
            data = data,
            mapping = aes_string(
                x = "sampleName",
                y = "nCells",
                fill = "interestingGroups"
            )
        ) +
            geom_bar(
                color = "black",
                stat = "identity"
            ) +
            labs(
                title = title,
                x = NULL,
                fill = paste(interestingGroups, collapse = ":\n")
            )

        # Color palette
        if (!is.null(fill)) {
            p <- p + fill
        }

        # Labels
        if (nrow(data) <= 16L) {
            p <- p + bcbio_geom_label(
                data = data,
                mapping = aes_string(label = "nCells"),
                # Align the label just under the top of the bar
                vjust = 1.25
            )
        }

        # Facets
        facets <- NULL
        if (.isAggregate(data)) {
            facets <- c(facets, "aggregate")
        }
        if (is.character(facets)) {
            p <- p + facet_wrap(facets = facets, scales = "free")
        }

        p
    }
)



#' @rdname seurat-SingleCellExperiment
#' @export
setMethod(
    "plotCellCounts",
    signature("seurat"),
    getMethod("plotCellCounts", "SingleCellExperiment")
)
