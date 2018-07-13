#' Quality Control Plots
#'
#' Utility function that loops our standard quality control plots, for easy
#' visualization.
#'
#' @name plotQC
#' @family Quality Control Functions
#' @author Michael Steinbaugh
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
#' plotQC(seurat_small)
NULL



#' @importFrom bcbioBase plotQC
#' @export
bcbioBase::plotQC



# Methods ======================================================================
#' @rdname plotQC
#' @export
setMethod(
    "plotQC",
    signature("SingleCellExperiment"),
    function(
        object,
        interestingGroups,
        geom = c("violin", "ridgeline", "ecdf", "histogram", "boxplot"),
        headerLevel = 2L,
        legend = FALSE,
        return = c("grid", "list", "markdown")
    ) {
        if (missing(interestingGroups)) {
            interestingGroups <- basejump::interestingGroups(object)
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
            "Cell Counts" = plotCellCounts,
            "Reads per Cell" = plotReadsPerCell,
            "UMIs per Cell" = plotUMIsPerCell,
            "Genes per Cell" = plotGenesPerCell,
            "UMIs vs. Genes" = plotUMIsVsGenes,
            "Novelty" = plotNovelty,
            "Mito Ratio" = plotMitoRatio,
            "Zeros vs. Depth" = plotZerosVsDepth
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
            plot_grid(plotlist = plotlist)
        } else if (return == "markdown") {
            markdownHeader(
                text = "Filtered quality control metrics",
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
