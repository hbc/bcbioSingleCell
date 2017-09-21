#' Plot Features of a Data Set
#'
#' @rdname plotFeatures
#' @name plotFeatures
#' @family Clustering Utilities
#' @author Michael Steinbaugh
#'
#' @param features Character vector of parameters supported by
#'   [Seurat::FetchData()] (e.g. nUMI, mitoRatio, PC).
#' @param headerLevel Include a Markdown header for each gene.
#'
#' @return No value, only graphical output.
NULL



# Constructors ====
.plotFeatureSeurat <- function(object, feature) {
    if (!is_string(feature)) {
        stop("feature must be a string", call. = FALSE)
    }

    # tSNE color plot
    # Dark theme enables greater contrast for marker visualization.
    # Otherwise use `rev(viridis(2))` for the colors. This will define
    # yellow as low and purple as high.
    tsne <- suppressMessages(FeaturePlot(
        object,
        features.plot = feature,
        # Use viridis for better contrast than default colors
        # (1) low: purple; (2) high: yellow
        do.return = TRUE,
        no.legend = FALSE) %>%
            # FeaturePlot outputs a named list - select the first item
            .[[1]] +
            # Use this for a better gradient
            scale_color_viridis(option = "inferno") +
            DarkTheme()
    )

    # Violin plot
    violin <- VlnPlot(
        object,
        features.plot = feature,
        do.return = TRUE,
        return.plotlist = TRUE) %>%
        .[[1]] %>%
        .[["data"]] %>%
        # Remove the low expression features
        dplyr::filter(feature > 0.1) %>%
        ggplot(aes_(x = ~ident,
                    y = ~feature,
                    fill = ~ident)) +
        geom_violin(
            color = NA,
            scale = "width",
            adjust = 1,
            trim = TRUE) +
        scale_fill_viridis(discrete = TRUE) +
        theme(legend.position = "none")

    # Ridgeline (joy) plot
    ridges <- JoyPlot(
        object,
        features.plot = feature,
        cols.use = viridis(length(levels(object@ident))),
        do.return = TRUE,
        return.plotlist = TRUE) %>%
        .[[1]] %>%
        .[["data"]] %>%
        # Remove the low expression features
        dplyr::filter(feature > 0.1) %>%
        ggplot(aes_(x = ~feature,
                    y = ~ident,
                    fill = ~ident)) +
        geom_density_ridges(color = NA, scale = 2) +
        scale_fill_viridis(discrete = TRUE) +
        theme(legend.position = "none")

    ggdraw() +
        # Coordinates are relative to lower left corner
        draw_plot(
            tsne,
            x = 0L, y = 0.3, width = 1, height = 0.7) +
        draw_plot(
            violin,
            x = 0, y = 0, width = 0.5, height = 0.3) +
        suppressMessages(draw_plot(
            ridges,
            x = 0.5, y = 0, width = 0.5, height = 0.3))
}



# Methods ====
#' @rdname plotFeatures
#' @export
setMethod("plotFeatures", "seurat", function(
    object,
    features,
    headerLevel = 2L) {
    lapply(seq_along(features), function(a) {
        mdHeader(features[[a]], level = headerLevel, asis = TRUE)
        .plotFeatureSeurat(object, features[[a]])
    }) %>%
        invisible
})
