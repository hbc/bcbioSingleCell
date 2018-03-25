#' Plot Mitochondrial Transcript Abundance
#'
#' @name plotMitoRatio
#' @family Quality Control Functions
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @inheritParams general
#'
#' @return `ggplot`.
#'
#' @examples
#' load(system.file("extdata/bcb_small.rda", package = "bcbioSingleCell"))
#' load(system.file("extdata/seurat_small.rda", package = "bcbioSingleCell"))
#'
#' # bcbioSingleCell ====
#' plotMitoRatio(bcb_small)
#'
#' # seurat ====
#' plotMitoRatio(seurat_small)
NULL



# Constructors =================================================================
.plotMitoRatio <- function(
    object,
    geom = c("violin", "boxplot", "histogram", "ridgeline"),
    interestingGroups,
    max,
    fill = scale_fill_viridis(discrete = TRUE)
) {
    geom <- match.arg(geom)
    .plotQCMetric(
        object = object,
        metricCol = "mitoRatio",
        geom = geom,
        interestingGroups = interestingGroups,
        max = max,
        trans = "identity",
        fill = fill
    )
}



# Methods ======================================================================
#' @rdname plotMitoRatio
#' @export
setMethod(
    "plotMitoRatio",
    signature("bcbioSingleCell"),
    .plotMitoRatio
)



#' @rdname plotMitoRatio
#' @export
setMethod(
    "plotMitoRatio",
    signature("seurat"),
    .plotMitoRatio
)
