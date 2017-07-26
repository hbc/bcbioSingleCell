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
    p <- ggplot(metrics,
           aes_(x = ~sample_name,
                y = ~mito_ratio * 100L,
                fill = as.name(interesting_group))) +
        labs(x = "sample",
             y = "% mito counts") +
        geom_boxplot() +
        geom_hline(alpha = qc_line_alpha,
                   color = qc_pass_color,
                   size = qc_line_size,
                   yintercept = max * 100L) +
        geom_label(data = median_mito_ratio,
                   aes_(label = ~round(mito_ratio * 100L,
                                       digits = 2L)),
                   alpha = qc_label_alpha,
                   label.padding = unit(0.1, "lines"),
                   show.legend = FALSE) +
        scale_y_sqrt() +
        theme(axis.text.x = element_text(angle = 90L, hjust = 1L))
    if (isTRUE(metadata(object)[["multiplexed_fastq"]])) {
        p <- p + facet_wrap(~file_name)
    }
    p
}



#' @rdname plot_mito_ratio
#' @usage NULL
.plot_mito_ratio_histogram <- function(object, max) {
    metrics <- metrics(object)
    p <- ggplot(metrics,
        aes_(x = ~mito_ratio * 100L,
             fill = ~sample_name)) +
        labs(x = "% mito counts") +
        geom_histogram(bins = bins) +
        geom_vline(alpha = qc_line_alpha,
                   color = qc_pass_color,
                   size = qc_line_size,
                   xintercept = max * 100L) +
        scale_x_sqrt() +
        scale_y_sqrt()
    if (isTRUE(metadata(object)[["multiplexed_fastq"]])) {
        p <- p + facet_wrap(~file_name)
    }
    p
}



#' @rdname plot_mito_ratio
#' @usage NULL
.plot_mito_ratio_scatterplot <- function(object) {
    metrics <- metrics(object)
    p <- ggplot(metrics,
        aes_(x = ~coding_counts,
             y = ~mito_counts,
             color = ~sample_name)) +
        labs(x = "mito counts",
             y = "coding counts") +
        geom_point(size = 1L) +
        scale_x_sqrt() +
        scale_y_sqrt() +
        theme(axis.text.x = element_text(angle = 90L, hjust = 1L))
    if (isTRUE(metadata(object)[["multiplexed_fastq"]])) {
        p <- p + facet_wrap(~file_name)
    }
    p
}



#' @rdname plot_mito_ratio
#' @usage NULL
.plot_mito_ratio <- function(object, max) {
    plot_grid(.plot_mito_ratio_scatterplot(object) +
                  theme(legend.position = "none"),
              .plot_mito_ratio_histogram(object, max) +
                  theme(legend.position = "none"),
              .plot_mito_ratio_boxplot(object, max) +
                  theme(legend.position = "bottom"),
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
