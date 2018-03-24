#' Plot Mitochondrial Counts vs. Coding Counts
#'
#' @name plotMitoVsCoding
#' @family Quality Control Functions
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @inheritParams general
#'
#' @examples
#' load(system.file("extdata/bcb_small.rda", package = "bcbioSingleCell"))
#' load(system.file("extdata/seurat_small.rda", package = "bcbioSingleCell"))
#'
#' # bcbioSingleCell ====
#' plotMitoVsCoding(bcb_small)
#'
#' # seurat ====
#' plotMitoVsCoding(seurat_small)
NULL



# Constructors =================================================================
.plotMitoVsCoding <- function(
    object,
    interestingGroups,
    color = scale_color_viridis(discrete = TRUE)
) {
    .plotQCScatterplot(
        object = object,
        xCol = "nCoding",
        yCol = "nMito",
        xTrans = "log2",
        yTrans = "log2"
    )
}



# Methods ======================================================================
#' @rdname plotMitoVsCoding
#' @export
setMethod(
    "plotMitoVsCoding",
    signature("bcbioSingleCell"),
    .plotMitoVsCoding
)



#' @rdname plotMitoVsCoding
#' @export
setMethod(
    "plotMitoVsCoding",
    signature("seurat"),
    .plotMitoVsCoding
)
