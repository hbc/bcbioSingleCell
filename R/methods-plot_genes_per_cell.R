#' Plot Genes per Cell
#'
#' @rdname plot_genes_per_cell
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @param min Recommended minimum value cutoff.
#' @param max Recommended maximum value cutoff.
#'
#' @return [ggplot].



#' @rdname plot_genes_per_cell
#' @usage NULL
.plot_genes_per_cell_boxplot <- function(object, min, max) {
    metrics <- metrics(object)
    median_genes <- aggregate(genes_detected ~ sample_id, metrics, median) %>%
        left_join(sample_metadata(object), by = "sample_id")
    interesting_group <- interesting_groups(object)[[1L]]
    plot <- ggplot(
        metrics,
        aes_(x = ~sample_name,
             y = ~genes_detected,
             fill = as.name(interesting_group))) +
        labs(x = "sample",
             y = "genes per cell") +
        facet_wrap(~file_name) +
        geom_boxplot() +
        geom_hline(alpha = qc_line_alpha,
                   color = qc_pass_color,
                   size = qc_line_size,
                   yintercept = min) +
        geom_label(data = median_genes,
                   aes_(label = ~round(genes_detected)),
                   alpha = qc_label_alpha,
                   label.padding = unit(0.1, "lines"),
                   show.legend = FALSE) +
        scale_y_sqrt() +
        theme(axis.text.x = element_text(angle = 90L, hjust = 1L))
    if (!is.null(max)) {
        plot <- plot +
            geom_hline(alpha = qc_line_alpha,
                       color = qc_pass_color,
                       size = qc_line_size,
                       yintercept = max)
    }
    plot
}



#' @rdname plot_genes_per_cell
#' @usage NULL
.plot_genes_per_cell_histogram <- function(object, min, max) {
    metrics <- metrics(object)
    plot <- ggplot(
        metrics,
        aes_(x = ~genes_detected,
             fill = ~sample_name)) +
        labs(x = "genes per cell") +
        facet_wrap(~file_name) +
        geom_histogram(bins = bins) +
        geom_vline(alpha = qc_line_alpha,
                   color = qc_pass_color,
                   size = qc_line_size,
                   xintercept = min) +
        scale_x_sqrt() +
        scale_y_sqrt() +
        theme(axis.text.x = element_text(angle = 90L, hjust = 1L))
    if (!is.null(max)) {
        plot <- plot +
            geom_vline(alpha = qc_line_alpha,
                       color = qc_pass_color,
                       size = qc_line_size,
                       xintercept = max)
    }
    plot
}



#' @rdname plot_genes_per_cell
#' @usage NULL
.plot_genes_per_cell <- function(object, min, max) {
    plot_grid(.plot_genes_per_cell_histogram(object, min, max),
              .plot_genes_per_cell_boxplot(object, min, max),
              labels = "auto",
              nrow = 2L)
}



#' @rdname plot_genes_per_cell
#' @export
setMethod(
    "plot_genes_per_cell",
    "bcbioSCDataSet",
    function(object, min = 500L, max = NULL) {
        .plot_genes_per_cell(object, min, max)
    })



#' @rdname plot_genes_per_cell
#' @export
setMethod(
    "plot_genes_per_cell",
    "bcbioSCSubset",
    function(object) {
        min <- object %>%
            metadata %>%
            .[["filtering_criteria"]] %>%
            .[["min_genes"]]
        max <- object %>%
            metadata %>%
            .[["filtering_criteria"]] %>%
            .[["max_genes"]]
        .plot_genes_per_cell(object, min, max)
    })
