# FIXME Round median labels -- pass digits here?

#' Plot Genes per Cell
#'
#' @name plotGenesPerCell
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
#' plotGenesPerCell(bcb_small)
#'
#' # seurat ====
#' plotGenesPerCell(seurat_small)
NULL



# Constructors =================================================================
.plotGenesPerCell <- function(
    object,
    geom = c("violin", "boxplot", "histogram", "ridgeline"),
    min = 0L,
    max = Inf,
    interestingGroups,
    flip = TRUE,
    fill = scale_fill_viridis(discrete = TRUE)
) {
    if (missing(interestingGroups)) {
        interestingGroups <- bcbioBase::interestingGroups(object)
    }
    if (missing(min)) {
        min <- metadata(object)[["filterParams"]][["minGenes"]]
    }
    if (missing(max)) {
        max <- metadata(object)[["filterParams"]][["maxGenes"]]
    }
    .plotQCGeom(
        metrics = metrics(object, interestingGroups = interestingGroups),
        metricCol = "nGene",
        geom = geom,
        min = min,
        max = max,
        interestingGroups = interestingGroups,
        flip = flip,
        fill = fill
    )
}



# Methods ======================================================================
#' @rdname plotGenesPerCell
#' @export
setMethod(
    "plotGenesPerCell",
    signature("bcbioSingleCell"),
    .plotGenesPerCell
)



#' @rdname plotGenesPerCell
#' @export
setMethod(
    "plotGenesPerCell",
    signature("seurat"),
    .plotGenesPerCell
)
