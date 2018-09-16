#' Plot Percentage of Zeros vs. Library Depth
#'
#' This function helps us visualize the dropout rate.
#'
#' @name plotZerosVsDepth
#' @family Quality Control Functions
#' @author Rory Kirchner, Michael Steinbaugh
#' @importFrom basejump plotZerosVsDepth
#' @export
#'
#' @inheritParams general
#'
#' @return `ggplot`.
#'
#' @examples
#' plotZerosVsDepth(indrops_small)
NULL



#' @rdname plotZerosVsDepth
#' @export
setMethod(
    "plotZerosVsDepth",
    signature("SingleCellExperiment"),
    function(
        object,
        interestingGroups = NULL,
        color = getOption("bcbio.discrete.color", NULL),
        title = "zeros vs. depth"
    ) {
        validObject(object)
        interestingGroups <- matchInterestingGroups(
            object = object,
            interestingGroups = interestingGroups
        )
        interestingGroups(object) <- interestingGroups
        assertIsColorScaleDiscreteOrNULL(color)
        assertIsAStringOrNULL(title)

        data <- zerosVsDepth(object)

        p <- ggplot(
            data = as(data, "tbl_df"),
            mapping = aes(
                x = !!sym("depth"),
                y = !!sym("dropout"),
                color = !!sym("interestingGroups")
            )
        ) +
            geom_point(size = 0.8, alpha = 0.8) +
            expand_limits(y = c(0L, 1L)) +
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



#' @rdname plotZerosVsDepth
#' @export
setMethod(
    "zerosVsDepth",
    signature("SingleCellExperiment"),
    function(object) {
        data <- zerosVsDepth(counts(object))
        # Stash the rownames.
        data[["rowname"]] <- rownames(data)
        # Add sample ID column.
        data[["sampleID"]] <- cell2sample(object)
        # Join the sample data.
        sampleData <- sampleData(object)
        if (is.null(sampleData)) {
            sampleData <- unknownSampleData
        }
        sampleData[["sampleID"]] <- rownames(sampleData)
        data <- merge(
            x = data,
            y = sampleData,
            by = "sampleID",
            all.x = TRUE
        )
        rownames(data) <- data[["rowname"]]
        data[["rowname"]] <- NULL
        # Ensure we're returning in the correct order.
        data <- data[colnames(object), , drop = FALSE]
        data
    }
)
