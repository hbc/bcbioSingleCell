#' @name plotCellCounts
#' @author Michael Steinbaugh, Rory Kirchner
#' @inherit bioverbs::plotCellCounts
#' @inheritParams basejump::params
#' @examples
#' data(indrops)
#' plotCellCounts(indrops)
NULL



#' @importFrom bioverbs plotCellCounts
#' @aliases NULL
#' @export
bioverbs::plotCellCounts



plotCellCounts.SingleCellExperiment <-  # nolint
    function(
        object,
        interestingGroups = NULL,
        fill,
        title = "cell counts"
    ) {
        validObject(object)
        assert(
            isGGScale(fill, scale = "discrete", aes = "fill") || is.null(fill),
            isString(title) || is.null(title)
        )
        interestingGroups(object) <-
            matchInterestingGroups(object, interestingGroups)

        metrics <- metrics(object)
        sampleData <- sampleData(object)

        # Remove user-defined `nCells` column, if present.
        metrics[["nCells"]] <- NULL
        sampleData[["nCells"]] <- NULL

        sampleData <- sampleData %>%
            as_tibble(rownames = "sampleID") %>%
            mutate_all(as.factor)

        data <- metrics %>%
            group_by(!!sym("sampleID")) %>%
            summarise(nCells = n()) %>%
            left_join(sampleData, by = "sampleID")

        p <- ggplot(
            data = data,
            mapping = aes(
                x = !!sym("sampleName"),
                y = !!sym("nCells"),
                fill = !!sym("interestingGroups")
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
            p <- p + basejump_geom_label(
                data = data,
                mapping = aes(label = !!sym("nCells")),
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
            p <- p + facet_wrap(facets = syms(facets), scales = "free")
        }

        p
    }

formals(plotCellCounts.SingleCellExperiment)[["fill"]] <-
    formalsList[["fill.discrete"]]


#' @rdname plotCellCounts
#' @export
setMethod(
    f = "plotCellCounts",
    signature = signature("SingleCellExperiment"),
    definition = plotCellCounts.SingleCellExperiment
)
