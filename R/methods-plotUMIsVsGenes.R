#' Plot UMI and Gene Correlation
#'
#' @name plotUMIsVsGenes
#' @family Quality Control Metrics
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @inheritParams general
#'
#' @examples
#' load(system.file("extdata/bcb_small.rda", package = "bcbioSingleCell"))
#' load(system.file("extdata/seurat_small.rda", package = "bcbioSingleCell"))
#'
#' # bcbioSingleCell ====
#' plotUMIsVsGenes(bcb_small)
#'
#' # seurat ====
#' plotUMIsVsGenes(pbmc_small)
#' plotUMIsVsGenes(seurat_small)
NULL



# Constructors =================================================================
.plotUMIsVsGenes <- function(
    object,
    interestingGroups,
    color = scale_color_viridis(discrete = TRUE)
) {
    .plotQCScatterplot(
        object = object,
        xCol = "nUMI",
        yCol = "nGene",
        xTrans = "log2",
        yTrans = "log2"
    )
}



# Methods ======================================================================
#' @rdname plotUMIsVsGenes
#' @export
setMethod(
    "plotUMIsVsGenes",
    signature("bcbioSingleCell"),
    .plotUMIsVsGenes
)



#' @rdname plotUMIsVsGenes
#' @export
setMethod(
    "plotUMIsVsGenes",
    signature("seurat"),
    .plotUMIsVsGenes
)
