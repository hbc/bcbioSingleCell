#' Plot novelty score
#'
#' "Novelty" refers to log10 genes detected per count.
#'
#' @rdname plot_novelty
#' @inheritParams plot_cell_counts



#' @rdname plot_novelty
#' @usage NULL
.plot_novelty_boxplot <- function(object, min) {
    metrics <- metrics(object)
    median_novelty <-
        aggregate(log10_genes_per_umi ~ sample_id, metrics, median) %>%
        left_join(sample_metadata(object), by = "sample_id")
    interesting_group <- interesting_groups(object)[[1L]]
    ggplot(
        metrics,
        aes_(x = ~sample_name,
             y = ~log10_genes_per_umi,
             fill = as.name(interesting_group))) +
        labs(x = "sample",
             y = "log10 genes per umi") +
        facet_wrap(~file_name) +
        geom_boxplot() +
        geom_hline(alpha = 0.5,
                   color = warn_color,
                   size = 2L,
                   yintercept = min) +
        geom_label(data = median_novelty,
                   aes_(label = ~round(log10_genes_per_umi, digits = 2L)),
                   alpha = 0.75,
                   label.padding = unit(0.1, "lines"),
                   show.legend = FALSE) +
        scale_y_sqrt() +
        theme(axis.text.x = element_text(angle = 90L, hjust = 1L))
}



#' @rdname plot_novelty
#' @usage NULL
.plot_novelty_histogram <- function(object, min) {
    metrics <- metrics(object)
    ggplot(
        metrics,
        aes_(x = ~log10_genes_per_umi,
             fill = ~sample_name)) +
        labs(x = "log10 genes per umi") +
        facet_wrap(~file_name) +
        geom_histogram(bins = bins) +
        geom_vline(alpha = 0.5,
                   color = warn_color,
                   size = 2L,
                   xintercept = min) +
        scale_x_sqrt() +
        scale_y_sqrt() +
        theme(axis.text.x = element_text(angle = 90L, hjust = 1L))
}



#' @rdname plot_novelty
#' @usage NULL
.plot_novelty <- function(object, min) {
    plot_grid(.plot_novelty_histogram(object, min),
              .plot_novelty_boxplot(object, min),
              labels = "auto",
              nrow = 2L)
}



#' @rdname plot_novelty
#' @export
setMethod(
    "plot_novelty",
    "bcbioSCDataSet",
    function(object, min = 0.8) {
        .plot_novelty(object, min)
    })

#' @rdname plot_novelty
#' @export
setMethod(
    "plot_novelty",
    "SCSubset",
    function(object) {
        min <- object %>%
            metadata %>%
            .[["filtering_criteria"]] %>%
            .[["min_novelty"]]
        .plot_novelty(object, min)
    })
