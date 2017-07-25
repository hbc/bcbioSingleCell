#' Plot Mitochondrial Transcript Abundance
#'
#' @rdname plot_mito_ratio
#' @author Michael Steinbaugh, Rory Kirchner
#' @inherit plot_genes_per_cell



#' @rdname plot_mito_ratio
#' @usage NULL
.plot_mito_ratio_boxplot <- function(object, max) {
    metrics <- metrics(object)
    median_mito_ratio <-
        aggregate(mito_ratio ~ sample_id, metrics, median) %>%
        left_join(sample_metadata(object), by = "sample_id")
    interesting_group <- interesting_groups(object)[[1L]]
    ggplot(
        metrics,
        aes_(x = ~sample_name,
             y = ~mito_ratio,
             fill = as.name(interesting_group))) +
        labs(x = "sample",
             y = "relative mitochondrial abundance") +
        facet_wrap(~file_name) +
        geom_boxplot() +
        geom_hline(alpha = 0.5,
                   color = warn_color,
                   size = 2L,
                   yintercept = max) +
        geom_label(data = median_mito_ratio,
                   aes_(label = ~round(mito_ratio, digits = 2L)),
                   alpha = 0.75,
                   label.padding = unit(0.1, "lines"),
                   show.legend = FALSE) +
        scale_y_sqrt() +
        theme(axis.text.x = element_text(angle = 90L, hjust = 1L))
}



#' @rdname plot_mito_ratio
#' @usage NULL
.plot_mito_ratio_histogram <- function(object, max) {
    metrics <- metrics(object)
    ggplot(
        metrics,
        aes_(x = ~mito_ratio,
             fill = ~sample_name)) +
        labs(x = "relative mitochondrial abundance") +
        facet_wrap(~file_name) +
        geom_histogram(bins = bins) +
        geom_vline(alpha = 0.5,
                   color = warn_color,
                   size = 2L,
                   xintercept = max) +
        scale_x_sqrt() +
        scale_y_sqrt() +
        theme(axis.text.x = element_text(angle = 90L, hjust = 1L))
}



#' @rdname plot_mito_ratio
#' @usage NULL
.plot_mito_ratio_scatterplot <- function(object) {
    metrics <- metrics(object)
    ggplot(
        metrics,
        aes_(x = ~coding_counts,
             y = ~mito_counts,
             color = ~sample_name)) +
        labs(x = "mitochondrial counts",
             y = "coding counts") +
        facet_wrap(~file_name) +
        geom_point(size = 1L) +
        scale_x_sqrt() +
        scale_y_sqrt() +
        theme(axis.text.x = element_text(angle = 90L, hjust = 1L))
}



#' @rdname plot_mito_ratio
#' @usage NULL
.plot_mito_ratio <- function(object, max) {
    plot_grid(.plot_mito_ratio_histogram(object, max),
              .plot_mito_ratio_boxplot(object, max),
              .plot_mito_ratio_scatterplot(object),
              labels = "auto",
              nrow = 3L)
}



#' @rdname plot_mito_ratio
#' @export
setMethod(
    "plot_mito_ratio",
    "bcbioSCDataSet",
    function(object, max = 0.1) {
        .plot_mito_ratio(object, max)
    })

#' @rdname plot_mito_ratio
#' @export
setMethod(
    "plot_mito_ratio",
    "bcbioSCSubset",
    function(object) {
        max <- object %>%
            metadata %>%
            .[["filtering_criteria"]] %>%
            .[["max_mito_ratio"]]
        .plot_mito_ratio(object, max)
    })
