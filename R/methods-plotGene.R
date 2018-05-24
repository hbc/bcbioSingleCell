#' Plot Gene
#'
#' @name plotGene
#' @family Gene Expression Functions
#' @author Michael Steinbaugh
#'
#' @importFrom bcbioBase plotGene
#'
#' @inheritParams general
#' @param geom Plot type. Uses [match.arg()] to pick the type. Currently
#'   supports "`dot`" and "`violin`".
#'
#' @seealso
#' - [plotDot()], [Seurat::DotPlot()].
#' - [plotViolin()], [Seurat::VlnPlot()].
#' - [Seurat::RidgePlot()].
#'
#' @return `ggplot`.
#'
#' @examples
#' object <- seurat_small
#' genes <- head(rownames(object))
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
    signature("seurat"),
    function(
        object,
        genes,
        geom = c("dot", "violin"),
        legend = FALSE
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
