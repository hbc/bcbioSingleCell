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
#' Plot t-SNE Constructor Function
#' @param object [data.frame] returned from [fetchTSNEExpressionData()].
.plotTSNE <- function(
    object,
    interestingGroup = "ident",
    label = TRUE) {
    if (interestingGroup == "ident") {
        color <- "object.ident"
    } else {
        color <- interestingGroup
    }
    p <- ggplot(object,
                aes_string(x = "tSNE_1",
                           y = "tSNE_2",
                           color = color)) +
        # Alpha transparency helps distinguish superimposed points
        geom_point(alpha = 0.6) +
        DarkTheme()
    if (isTRUE(label)) {
        p <- p +
            geom_text(
                aes_string(
                    x = "centerx",
                    y = "centery",
                    label = "object.ident"),
                color = "white",
                size = 6,
                fontface = "bold")
    }
    p <- p +
        # Use the default color palette here -- more dynamic range than viridis
        guides(color = guide_legend(
            title.position = "left",
            nrow = 2,
            byrow = TRUE)) +
        theme(legend.justification = "center",
              legend.position = "bottom")
    if (interestingGroup == "ident") {
        # Fix the cluster identity label
        p <- p + labs(color = "cluster")
    }
    p
}



# Methods ====
#' @rdname plotTSNE
#' @export
setMethod("plotTSNE", "seurat", function(
    object,
    interestingGroup = "ident") {
    tsne <- fetchTSNEData(object)
    .plotTSNE(tsne, interestingGroup = interestingGroup)
})
