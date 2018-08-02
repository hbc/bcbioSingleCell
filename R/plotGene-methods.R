#' Plot Gene
#'
#' Visualize genes on a dot or violin plot.
#'
#' @name plotGene
#' @family Gene Expression Functions
#' @author Michael Steinbaugh
#'
#' @importFrom bcbioBase plotGene
#' @export
#'
#' @inheritParams general
#' @param geom `string`. Plot type. Uses [match.arg()] to pick the type.
#'   Currently supports "`dot`" and "`violin`".
#'
#' @seealso
#' - [plotDot()], [Seurat::DotPlot()].
#' - [plotViolin()], [Seurat::VlnPlot()].
#' - [Seurat::RidgePlot()].
#'
#' @return `ggplot`.
#'
#' @examples
#' object <- indrops_small
#' genes <- head(rownames(object))
#' glimpse(genes)
#'
#' # Dot
#' plotGene(object, genes = genes, geom = "dot")
#'
#' # Violin
#' plotGene(object, genes = genes, geom = "violin")
NULL



# Methods ======================================================================
#' @rdname plotGene
#' @export
setMethod(
    "plotGene",
    signature("SingleCellExperiment"),
    function(
        object,
        genes,
        geom = c("dot", "violin"),
        legend = getOption("bcbio.legend", TRUE)
    ) {
        geom <- match.arg(geom)
        if (geom == "dot") {
            fun <- plotDot
        } else if (geom == "violin") {
            fun <- plotViolin
        }
        fun(
            object = object,
            genes = genes,
            legend = legend
        )
    }
)
