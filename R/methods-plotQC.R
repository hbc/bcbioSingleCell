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
#' plotQC(bcb_small)
#'
#' # SingleCellExperiment ====
#' suppressWarnings(plotQC(cellranger_small))
#'
#' # seurat ====
#' suppressWarnings(plotQC(Seurat::pbmc_small))
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
        geom = c("violin", "boxplot", "histogram", "ridgeline"),
        headerLevel = 2L,
        legend = FALSE,
        return = c("grid", "list", "markdown")
    ) {
        if (missing(interestingGroups)) {
            interestingGroups <- bcbioBase::interestingGroups(object)
        }
        geom <- match.arg(geom)
        return <- match.arg(return)

        plotlist <- list(
            plotReadsPerCell = plotReadsPerCell(
                object,
                interestingGroups = interestingGroups
            ),
            plotCellCounts = plotCellCounts(
                object,
                interestingGroups = interestingGroups
            ),
            plotUMIsPerCell = plotUMIsPerCell(
                object,
                interestingGroups = interestingGroups,
                geom = geom
            ),
            plotGenesPerCell = plotGenesPerCell(
                object,
                interestingGroups = interestingGroups,
                geom = geom
            ),
            plotUMIsVsGenes = plotUMIsVsGenes(
                object,
                interestingGroups = interestingGroups
            ),
            plotMitoRatio = plotMitoRatio(
                object,
                interestingGroups = interestingGroups,
                geom = geom
            ),
            plotNovelty = plotNovelty(
                object,
                interestingGroups = interestingGroups,
                geom = geom
            )
        )

        # Remove any `NULL` plots. This is useful for nuking the
        # `plotReadsPerCell()` return on an object that doesn't contain raw cellular
        # barcode counts.
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

            # Reads per cell currently only supported for bcbio runs
            if (!is.null(plotlist[["plotReadsPerCell"]])) {
                markdownHeader(
                    "Reads per cell",
                    level = headerLevel + 1L,
                    asis = TRUE
                )
                show(plotlist[["plotReadsPerCell"]])
            }

            markdownHeader(
                "Cell counts",
                level = headerLevel + 1L,
                asis = TRUE
            )
            show(plotlist[["plotCellCounts"]])

            markdownHeader(
                "UMI counts per cell",
                level = headerLevel + 1L,
                asis = TRUE
            )
            show(plotlist[["plotUMIsPerCell"]])

            markdownHeader(
                "Genes detected",
                level = headerLevel + 1L,
                asis = TRUE
            )
            show(plotlist[["plotGenesPerCell"]])

            markdownHeader(
                "UMIs vs. genes",
                level = headerLevel + 1L,
                asis = TRUE
            )
            show(plotlist[["plotUMIsVsGenes"]])

            markdownHeader(
                "Mitochondrial counts ratio",
                level = headerLevel + 1L,
                asis = TRUE
            )
            show(plotlist[["plotMitoRatio"]])

            markdownHeader(
                "Novelty",
                level = headerLevel + 1L,
                asis = TRUE
            )
            show(plotlist[["plotNovelty"]])
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
