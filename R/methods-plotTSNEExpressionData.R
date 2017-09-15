#' Plot Fetched t-SNE Expression Data
#'
#' @rdname plotTSNEExpressionData
#' @name plotTSNEExpressionData
#' @family t-SNE Utilities
#' @author Rory Kirchner, Michael Steinbaugh
#'
#' @param object [data.frame] returned from [fetchTSNEExpressionData()].
#' @param colorpoints Color points by geometric mean or expression of
#'  individual gene.
#'
#' @return t-SNE [ggplot].
#'
#' @examples
#' \dontrun{
#' data(markers, seurat)
#' top <- topMarkers(markers)
#'
#' # Let's take the top markers specific to cluster 0, as an example
#' cluster0 <- top %>%
#'     filter(cluster == 0) %>%
#'     pull(symbol)
#'
#' # Fetch the t-SNE expression data for the desired gene symbols
#' dat <- fetchTSNEExpressionData(seurat, cluster0)
#' print(unique(dat$gene))
#'
#' # To make t-SNE plot colored by geometric mean of topGenes
#' plotTSNEExpressionData(dat, colorpoints = "geomean")
#'
#' # To make faceted t-SNE plot of each gene (looks good at up to 6 genes)
#' plotTSNEExpressionData(dat, colorpoints = "expression") +
#'     facet_wrap(~gene)
#' }
setMethod("plotTSNEExpressionData", "grouped_df", function(
    object,
    colorpoints = "geomean") {
    if (!colorpoints %in% c("expression", "geomean")) {
        stop("colorpoints supports 'geomean' or 'expression'")
    }
    # Check to make sure only a subset of markers is passed in
    if (nrow(object) > 30) {
        stop("'plotTSNEExpressionData()' supports up to 30 marker genes")
    }
    object %>%
        ggplot(
            aes_string(
                x = "tSNE_1",
                y = "tSNE_2",
                color = colorpoints)) +
        geom_text(
            aes_string(
                x = "tSNE_1",
                y = "tSNE_2",
                label = "object.ident",
                color = colorpoints),
            size = 2) +
        scale_color_viridis() +
        guides(color = FALSE) +
        geom_text(
            aes_string(
                x = "centerx",
                y = "centery",
                label = "object.ident"),
            color = "pink",
            size = 3,
            fontface = "bold") +
        xlab("") +
        ylab("") +
        DarkTheme() +
        theme(strip.background = element_rect(fill = "black"),
              strip.text = element_text(face = "bold"),
              text = element_text(size = 8))
})
