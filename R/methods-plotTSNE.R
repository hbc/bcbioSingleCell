#' Plot t-SNE
#'
#' Generate a t-distributed stochastic neighbor embedding (t-SNE) plot.
#'
#' @rdname plotTSNE
#' @name plotTSNE
#' @family t-SNE Utilities
#' @author Rory Kirchner, Michael Steinbaugh
#'
#' @inheritParams AllGenerics
#' @param interestingGroup Interesting group to use for plot colors.
#' @param label Label the clusters on the plot.
#'
#' @return [ggplot].
#'
#' @examples
#' \dontrun{
#' plotTSNE(seurat)
#' }
NULL



# Constructors ====
#' Plot t-SNE Constructor
#'
#' @author Rory Kirchner, Michael Steinbaugh
#'
#' @param object [data.frame] returned from [fetchTSNEExpressionData()].
#'
#' @return [ggplot].
#' @noRd
.plotDim <- function(
    object,
    axes,
    interestingGroup = "ident",
    label = TRUE) {
    if (interestingGroup == "ident") {
        # Seurat stores the ident from `FetchData()` as `object.ident`
        color <- "ident"
    } else {
        color <- interestingGroup
    }
    p <- ggplot(
        object,
        mapping = aes_string(
            x = axes[["x"]],
            y = axes[["y"]],
            color = color)
    ) +
        # Alpha transparency helps distinguish superimposed points
        geom_point(alpha = 0.7) +
        darkTheme() +
        guides(color = guide_legend(
            title.position = "left",
            byrow = TRUE))
    if (isTRUE(label)) {
        p <- p +
            geom_text(
                mapping = aes_string(
                    x = "centerX",
                    y = "centerY",
                    label = "ident"),
                color = "white",
                size = 6,
                fontface = "bold")
    }
    if (interestingGroup == "ident") {
        # Fix the cluster identity label
        p <- p +
            labs(color = "cluster")
    }
    p
}



# Methods ====
#' @rdname plotTSNE
#' @export
setMethod("plotTSNE", "seurat", function(
    object,
    interestingGroup = "ident",
    label = TRUE) {
    tsne <- fetchTSNEData(object)
    .plotDim(
        tsne,
        axes = c(x = "tSNE1", y = "tSNE2"),
        interestingGroup = interestingGroup,
        label = label)
})
