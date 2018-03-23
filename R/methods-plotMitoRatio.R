#' Plot Mitochondrial Transcript Abundance
#'
#' @name plotMitoRatio
#' @family Quality Control Metrics
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
    max = Inf,
    interestingGroups,
    flip = TRUE,
    fill = scale_fill_viridis(discrete = TRUE)
) {
    if (missing(interestingGroups)) {
        interestingGroups <- bcbioBase::interestingGroups(object)
    }
    if (missing(max)) {
        max <- metadata(object)[["filterParams"]][["maxMitoRatio"]]
    }
    .plotQCGeom(
        metrics = metrics(object, interestingGroups = interestingGroups),
        metricCol = "mitoRatio",
        geom = geom,
        max = max,
        interestingGroups = interestingGroups,
        flip = flip,
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
