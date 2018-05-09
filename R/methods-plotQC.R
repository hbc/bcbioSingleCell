#' Quality Control Plots
#'
#' Utility function that loops our standard quality control plots, for easy
#' visualization.
#'
#' @name plotQC
#' @family Quality Control Functions
#' @author Michael Steinbaugh
#'
#' @importFrom bcbioBase plotQC
#'
#' @inheritParams general
#'
#' @return
#'   - `grid`: [cowplot::plot_grid()] graphical output.
#'   - `list`: `list` containing `ggplot` objects.
#'   - `markdown`: R Markdown report, with reports separated by headers.
#'
#' @examples
#' # bcbioSingleCell ====
#' plotQC(indrops_small)
#'
#' # SingleCellExperiment ====
#' plotQC(cellranger_small)
#'
#' # seurat ====
#' plotQC(Seurat::pbmc_small)
NULL



# Methods ======================================================================
#' @rdname plotQC
#' @export
setMethod(
    "plotQC",
    signature("SingleCellExperiment"),
    function(
        object,
        interestingGroups,
        geom = c("ecdf", "histogram", "ridgeline", "violin", "boxplot"),
        headerLevel = 2L,
        legend = FALSE,
        return = c("grid", "list", "markdown")
    ) {
        if (missing(interestingGroups)) {
            interestingGroups <- bcbioBase::interestingGroups(object)
        } else {
            interestingGroups(object) <- interestingGroups
        }
        geom <- match.arg(geom)
        return <- match.arg(return)

        plotCellCounts <- plotCellCounts(object)
        plotReadsPerCell <- NULL
        plotUMIsPerCell <- plotUMIsPerCell(object, geom = geom)
        plotGenesPerCell <- plotGenesPerCell(object, geom = geom)
        plotUMIsVsGenes <- plotUMIsVsGenes(object)
        plotNovelty <- plotNovelty(object, geom = geom)
        plotMitoRatio <- plotMitoRatio(object, geom = geom)
        plotZerosVsDepth <- plotZerosVsDepth(object)

        if (is(object, "bcbioSingleCell")) {
            # Don't show cell counts for unfiltered bcbio datasets
            if (!length(metadata(object)[["filterCells"]])) {
                plotCellCounts <- NULL
            }
            # Raw read counts are only stashed in bcbioSingleCell objects
            plotReadsPerCell <- plotReadsPerCell(object, geom = geom)
        }

        plotlist <- list(
            "plotCellCounts" = plotCellCounts,
            "plotReadsPerCell" = plotReadsPerCell,
            "plotUMIsPerCell" = plotUMIsPerCell,
            "plotGenesPerCell" = plotGenesPerCell,
            "plotUMIsVsGenes" = plotUMIsVsGenes,
            "plotNovelty" = plotNovelty,
            "plotMitoRatio" = plotMitoRatio,
            "plotZerosVsDepth" = plotZerosVsDepth
        )

        # Remove any `NULL` plots. This is useful for nuking the
        # `plotReadsPerCell()` return on an object that doesn't contain raw
        # cellular barcode counts.
        plotlist <- Filter(Negate(is.null), plotlist)

        # Hide the legends, if desired
        if (identical(legend, FALSE)) {
            .hideLegend <- function(gg) {
                gg + theme(legend.position = "none")
            }
            plotlist <- lapply(plotlist, .hideLegend)
        }

        # Grid return mode
        if (return == "list") {
            plotlist
        } else if (return == "grid") {
            plot_grid(plotlist = plotlist, labels = "AUTO")
        } else if (return == "markdown") {
            markdownHeader(
                "Filtered quality control metrics",
                level = headerLevel,
                tabset = TRUE,
                asis = TRUE
            )
            markdownPlotlist(
                plotlist = plotlist,
                headerLevel = headerLevel + 1L
            )
        }
    }
)



#' @rdname plotQC
#' @export
setMethod(
    "plotQC",
    signature("seurat"),
    getMethod("plotQC", "SingleCellExperiment")
)
