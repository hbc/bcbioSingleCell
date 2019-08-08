#' @name plotQC
#' @author Michael Steinbaugh
#' @include globals.R
#' @inherit acidplots::plotQC
#' @note Updated 2019-07-29.
#'
#' @inheritParams acidroxygen::params
#' @param headerLevel `integer(1)` (`1`-`7`).
#'   R Markdown header level.
#' @param legend `logical(1)`.
#'   Include plot legend.
#' @param ... Additional arguments.
#'
#' @return
#'   - `grid`: Graphical output.
#'   - `list`: `list` containing `ggplot` objects.
#'   - `markdown`: R Markdown report, with reports separated by headers.
#'
#' @examples
#' data(indrops)
#' plotQC(indrops)
NULL



#' @rdname plotQC
#' @name plotQC
#' @importFrom bioverbs plotQC
#' @usage plotQC(object, ...)
#' @export
NULL



## Updated 2019-07-24.
`plotQC,bcbioSingleCell` <-  # nolint
    function(
        object,
        interestingGroups = NULL,
        geom,
        headerLevel = 2L,
        legend,
        return = c("grid", "list", "markdown")
    ) {
        validObject(object)
        assert(isHeaderLevel(headerLevel))
        interestingGroups(object) <-
            matchInterestingGroups(object, interestingGroups)
        geom <- match.arg(geom)
        return <- match.arg(return)

        ## Don't show cell counts for unfiltered datasets.
        if (!is.null(metadata(object)[["filterCells"]])) {
            plotCellCounts <- plotCellCounts(object)
            plotZerosVsDepth <- NULL
        } else {
            plotCellCounts <- NULL
            plotZerosVsDepth <- plotZerosVsDepth(object)
        }

        plotUMIsPerCell <- plotUMIsPerCell(object, geom = geom)
        plotGenesPerCell <- plotGenesPerCell(object, geom = geom)
        plotUMIsVsGenes <- plotUMIsVsGenes(object)
        plotNovelty <- plotNovelty(object, geom = geom)
        plotMitoRatio <- plotMitoRatio(object, geom = geom)

        list <- list(
            "Cell counts" = plotCellCounts,
            "UMIs per cell" = plotUMIsPerCell,
            "Genes per cell" = plotGenesPerCell,
            "UMIs vs. genes" = plotUMIsVsGenes,
            "Novelty" = plotNovelty,
            "Mito ratio" = plotMitoRatio,
            "Zeros vs. depth" = plotZerosVsDepth
        )

        ## Remove any `NULL` plots. This is useful for nuking the
        ## `plotReadsPerCell` return on an object that doesn't contain raw
        ## cellular barcode counts.
        list <- Filter(f = Negate(is.null), x = list)

        ## Consistently show n plots.
        n <- 6L
        assert(hasLength(list, n = n))

        ## Hide the legends, if desired.
        if (identical(legend, FALSE)) {
            .hideLegend <- function(gg) {
                gg + theme(legend.position = "none")
            }
            list <- lapply(list, .hideLegend)
        }

        ## Return.
        if (return == "list") {
            names(list) <- camelCase(names(list))
            list
        } else if (return == "grid") {
            plot_grid(
                plotlist = list,
                ncol = n / 2L,
                nrow = 2L
            )
        } else if (return == "markdown") {
            markdownHeader(
                text = "Quality control metrics",
                level = headerLevel,
                tabset = TRUE,
                asis = TRUE
            )
            markdownPlots(list = list, headerLevel = headerLevel + 1L)
        }
    }

formals(`plotQC,bcbioSingleCell`)[["geom"]] <- geom
formals(`plotQC,bcbioSingleCell`)[["legend"]] <- formalsList[["legend"]]



#' @rdname plotQC
#' @export
setMethod(
    f = "plotQC",
    signature = signature("bcbioSingleCell"),
    definition = `plotQC,bcbioSingleCell`
)
