#' @name plotCellCounts
#' @author Michael Steinbaugh, Rory Kirchner
#' @inherit bioverbs::plotCellCounts
#' @note Updated 2019-07-29.
#'
#' @inheritParams acidroxygen::params
#' @param ... Additional arguments.
#'
#' @examples
#' data(indrops)
#' plotCellCounts(indrops)
NULL



#' @rdname plotCellCounts
#' @name plotCellCounts
#' @importFrom bioverbs plotCellCounts
#' @usage plotCellCounts(object, ...)
#' @export
NULL



## Updated 2019-07-27.
`plotCellCounts,bcbioSingleCell` <-  # nolint
    function(
        object,
        interestingGroups = NULL,
        fill,
        title = "Cell counts"
    ) {
        validObject(object)
        assert(
            isGGScale(fill, scale = "discrete", aes = "fill", nullOK = TRUE),
            isString(title, nullOK = TRUE)
        )
        interestingGroups(object) <-
            matchInterestingGroups(object, interestingGroups)

        metrics <- metrics(object)
        sampleData <- sampleData(object)
        metricCol <- "nCells"

        ## Remove user-defined column, if present.
        metrics[[metricCol]] <- NULL
        sampleData[[metricCol]] <- NULL

        sampleData <- sampleData %>%
            as_tibble(rownames = "sampleID") %>%
            mutate_all(as.factor)

        data <- metrics %>%
            group_by(!!sym("sampleID")) %>%
            summarise(!!sym(metricCol) := n()) %>%
            left_join(sampleData, by = "sampleID")

        p <- ggplot(
            data = data,
            mapping = aes(
                x = !!sym("sampleName"),
                y = !!sym(metricCol),
                fill = !!sym("interestingGroups")
            )
        ) +
            acid_geom_bar() +
            acid_scale_y_continuous_nopad() +
            labs(
                title = title,
                x = NULL,
                y = makeLabel(metricCol),
                fill = paste(interestingGroups, collapse = ":\n")
            )

        ## Color palette
        if (!is.null(fill)) {
            p <- p + fill
        }

        ## Labels
        if (nrow(data) <= 16L) {
            p <- p + acid_geom_label(
                data = data,
                mapping = aes(label = !!sym(metricCol)),
                ## Align the label just under the top of the bar.
                vjust = 1.25
            )
        }

        ## Facets
        facets <- NULL
        if (.isAggregate(data)) {
            facets <- c(facets, "aggregate")
        }
        if (is.character(facets)) {
            p <- p + facet_wrap(facets = syms(facets), scales = "free")
        }

        p
    }

formals(`plotCellCounts,bcbioSingleCell`)[["fill"]] <-
    formalsList[["fill.discrete"]]


#' @rdname plotCellCounts
#' @export
setMethod(
    f = "plotCellCounts",
    signature = signature("bcbioSingleCell"),
    definition = `plotCellCounts,bcbioSingleCell`
)
