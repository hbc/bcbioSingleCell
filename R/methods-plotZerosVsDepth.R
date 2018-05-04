#' Plot Percentage of Zeros vs. Library Depth
#'
#' This function helps us visualize the dropout rate.
#'
#' @name plotZerosVsDepth
#' @family Quality Control Functions
#' @author Rory Kirchner, Michael Steinbaugh
#'
#' @inheritParams general
#'
#' @return `ggplot`.
#'
#' @examples
#' # bcbioSingleCell ====
#' plotZerosVsDepth(indrops_small)
#'
#' # SingleCellExperiment ====
#' plotZerosVsDepth(cellranger_small)
#'
#' # seurat ====
#' plotZerosVsDepth(Seurat::pbmc_small)
NULL



# Methods ======================================================================
#' @rdname plotZerosVsDepth
#' @export
setMethod(
    "plotZerosVsDepth",
    signature("SingleCellExperiment"),
    function(
        object,
        interestingGroups,
        color = scale_color_hue(),
        title = "zeros vs. depth"
    ) {
        if (missing(interestingGroups)) {
            interestingGroups <- bcbioBase::interestingGroups(object)
        }
        assertIsColorScaleDiscreteOrNULL(color)
        assertIsAStringOrNULL(title)

        sampleData <- sampleData(
            object = object,
            interestingGroups = interestingGroups,
            return = "data.frame"
        )
        sampleData[["sampleID"]] <- as.factor(rownames(sampleData))

        counts <- assay(object)
        # Using a logical matrix is faster and more memory efficient
        present <- counts %>%
            # Ensure dgTMatrix gets coereced (e.g. Seurat::pbmc_small)
            as("dgCMatrix") %>%
            as("lgCMatrix")

        data <- tibble(
            "sampleID" = cell2sample(object),
            "dropout" = (nrow(present) - colSums(present)) / nrow(present),
            "depth" = colSums(counts)
        ) %>%
            left_join(sampleData, by = "sampleID")

        p <- ggplot(
            data = data,
            mapping = aes_string(
                x = "depth",
                y = "dropout",
                color = "interestingGroups"
            )
        ) +
            geom_point(size = 0.8, alpha = 0.8) +
            scale_x_continuous(trans = "log10") +
            labs(
                title = title,
                x = "library size (depth)",
                y = "dropout rate",
                color = paste(interestingGroups, collapse = ":\n")
            )

        if (is(color, "ScaleDiscrete")) {
            p <- p + color
        }

        # Wrap aggregated samples
        facets <- NULL
        if (.isAggregate(data)) {
            facets <- "aggregate"
        }
        if (is.character(facets)) {
            p <- p + facet_wrap(facets = facets, scales = "free")
        }

        p
    }
)



#' @rdname plotZerosVsDepth
#' @export
setMethod(
    "plotZerosVsDepth",
    signature("seurat"),
    getMethod("plotZerosVsDepth", "SingleCellExperiment")
)
