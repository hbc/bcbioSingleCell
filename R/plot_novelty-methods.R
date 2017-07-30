#' Plot Novelty Score
#'
#' @details "Novelty" refers to log10 genes detected per count.
#'
#' @rdname plot_novelty
#' @inherit plot_genes_per_cell



#' @rdname plot_novelty
#' @usage NULL
.plot_novelty_boxplot <- function(object, min) {
    metrics <- metrics(object)
    median_novelty <-
        aggregate(log10_genes_per_umi ~ sampleID, metrics, median) %>%
        left_join(sample_metadata(object), by = "sampleID")
    interesting_group <- interesting_groups(object)[[1L]]
    p <- ggplot(metrics,
        aes_(x = ~sample_name,
             y = ~log10_genes_per_umi,
             fill = as.name(interesting_group))) +
        labs(x = "sample",
             y = "log10 genes per umi") +
        geom_boxplot() +
        geom_hline(alpha = qc_line_alpha,
                   color = qc_pass_color,
                   size = qc_line_size,
                   yintercept = min) +
        geom_label(data = median_novelty,
                   aes_(label = ~round(log10_genes_per_umi, digits = 2L)),
                   alpha = qc_label_alpha,
                   label.padding = unit(0.1, "lines"),
                   show.legend = FALSE) +
        scale_y_sqrt() +
        theme(axis.text.x = element_text(angle = 90L, hjust = 1L))
    if (isTRUE(metadata(object)[["multiplexedFASTQ"]])) {
        p <- p + facet_wrap(~file_name)
    }
    p
}



#' @rdname plot_novelty
#' @usage NULL
.plot_novelty_histogram <- function(object, min) {
    metrics <- metrics(object)
    p <- ggplot(metrics,
        aes_(x = ~log10_genes_per_umi,
             fill = ~sample_name)) +
        labs(x = "log10 genes per umi") +
        geom_histogram(bins = bins) +
        geom_vline(alpha = qc_line_alpha,
                   color = qc_pass_color,
                   size = qc_line_size,
                   xintercept = min) +
        scale_x_sqrt() +
        scale_y_sqrt()
    if (isTRUE(metadata(object)[["multiplexedFASTQ"]])) {
        p <- p + facet_wrap(~file_name)
    }
    p
}



#' @rdname plot_novelty
#' @usage NULL
.plot_novelty <- function(object, min) {
    plot_grid(.plot_novelty_histogram(object, min) +
                  theme(legend.position = "none"),
              .plot_novelty_boxplot(object, min) +
                  theme(legend.position = "bottom"),
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
    "bcbioSCSubset",
    function(object) {
        min <- object %>%
            metadata %>%
            .[["filtering_criteria"]] %>%
            .[["min_novelty"]]
        .plot_novelty(object, min)
    })
