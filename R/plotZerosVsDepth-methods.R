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
#' plotZerosVsDepth(indrops_small)
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
        color = getOption("bcbio.discrete.color", NULL),
        title = "zeros vs. depth"
    ) {
        validObject(object)
        interestingGroups <- .prepareInterestingGroups(
            object = object,
            interestingGroups = interestingGroups
        )
        assertIsColorScaleDiscreteOrNULL(color)
        assertIsAStringOrNULL(title)

        counts <- counts(object)
        # Using a logical matrix is faster and more memory efficient
        present <- counts %>%
            # Ensure dgTMatrix gets coereced to dgCMatrix prior to logical
            as("dgCMatrix") %>%
            as("lgCMatrix")
        data <- tibble(
            sampleID = cell2sample(object),
            dropout = (nrow(present) - colSums(present)) / nrow(present),
            depth = colSums(counts)
        )

        sampleData <- sampleData(
            object = object,
            interestingGroups = interestingGroups
        )
        if (is.null(sampleData)) {
            sampleData <- unknownSampleData
        }
        sampleData <- as.data.frame(sampleData)
        sampleData[["sampleID"]] <- factor(
            x = rownames(sampleData),
            levels = levels(data[["sampleID"]])
        )

        data <- left_join(data, sampleData, by = "sampleID")

        p <- ggplot(
            data = data,
            mapping = aes(
                x = !!sym("depth"),
                y = !!sym("dropout"),
                color = !!sym("interestingGroups")
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
            p <- p + facet_wrap(facets = syms(facets), scales = "free")
        }

        p
    }
)
