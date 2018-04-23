#' Plot UMI and Gene Correlation
#'
#' @name plotUMIsVsGenes
#' @family Quality Control Metrics
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @inheritParams general
#'
#' @examples
#' # bcbioSingleCell ====
#' plotUMIsVsGenes(bcb_small)
#'
#' # SingleCellExperiment ====
#' plotUMIsVsGenes(cellranger_small)
#'
#' # seurat ====
#' plotUMIsVsGenes(Seurat::pbmc_small)
NULL



# Methods ======================================================================
#' @rdname plotUMIsVsGenes
#' @export
setMethod(
    "plotUMIsVsGenes",
    signature("SingleCellExperiment"),
    function(
        object,
        interestingGroups,
        color = scale_color_hue()
    ) {
        .plotQCScatterplot(
            object = object,
            xCol = "nUMI",
            yCol = "nGene",
            xTrans = "log2",
            yTrans = "log2"
        )
    }
)



#' @rdname plotUMIsVsGenes
#' @export
setMethod(
    "plotUMIsVsGenes",
    signature("seurat"),
    getMethod("plotUMIsVsGenes", "SingleCellExperiment")
)
