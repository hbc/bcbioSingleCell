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
#' # SingleCellExperiment ====
#' plotZerosVsDepth(bcb_small)
#' plotZerosVsDepth(cellranger_small)
NULL



# Methods ======================================================================
#' @rdname plotZerosVsDepth
#' @export
setMethod(
    "plotZerosVsDepth",
    signature("SingleCellExperiment"),
    function(
        object,
        color = scale_color_viridis(discrete = TRUE)
    ) {
        assertIsColorScaleDiscreteOrNULL(color)

        counts <- assay(object)
        sampleID <- cell2sample(object)
        sampleData <- sampleData(object, return = "data.frame")

        # Using a logical matrix is faster and more memory efficient
        present <- counts %>%
            # Ensure dgTMatrix gets coereced (e.g. pbmc_small)
            as("dgCMatrix") %>%
            as("lgCMatrix")

        data <- tibble(
            "sampleID" = sampleID,
            "dropout" = (nrow(present) - Matrix::colSums(present)) / nrow(present),
            "depth" = Matrix::colSums(counts)
        ) %>%
            left_join(sampleData, by = "sampleID")

        p <- ggplot(
            data = data,
            mapping = aes_string(
                x = "depth",
                y = "dropout",
                color = "sampleName"
            )
        ) +
            geom_point(size = 0.8, alpha = 0.8) +
            scale_x_continuous(trans = "log10") +
            labs(x = "library size (depth)", y = "dropout rate")

        if (is(color, "ScaleDiscrete")) {
            p <- p + color
        }

        # Wrap aggregated samples
        facets <- NULL
        if (.isAggregate(metrics)) {
            facets <- "sampleNameAggregate"
        }
        if (is.character(facets)) {
            p <- p + facet_wrap(facets = facets, scales = "free")
        }

        p
    }
)
