#' Plot UMIs per cell
#'
#' Plot the universal molecular identifiers (UMIs) per cell.
#'
#' @rdname plot_umis_per_cell
#' @author Michael Steinbaugh, Rory Kirchner
#' @inherit plot_genes_per_cell




#' @rdname plot_umis_per_cell
#' @usage NULL
.plot_umis_per_cell_boxplot <- function(object, min) {
    metrics <- metrics(object)
    median_umis <- aggregate(umi_counts ~ sample_id, metrics, median) %>%
        left_join(sample_metadata(object), by = "sample_id") %>%
        mutate(umi_counts = round(.data[["umi_counts"]]))
    interesting_group <- interesting_groups(object)[[1L]]
    ggplot(
        metrics,
        aes_(x = ~sample_name,
             y = ~umi_counts,
             fill = as.name(interesting_group))) +
        labs(title = "umi counts boxplot",
             x = "sample",
             y = "umis per cell") +
        facet_wrap(~file_name) +
        geom_boxplot() +
        geom_label(data = median_umis,
                   aes_(label = ~umi_counts),
                   alpha = 0.75,
                   label.padding = unit(0.1, "lines"),
                   show.legend = FALSE) +
        geom_hline(alpha = 0.5,
                   color = warn_color,
                   size = 2L,
                   yintercept = min) +
        scale_y_log10() +
        theme(axis.text.x = element_text(angle = 90L, hjust = 1L))
}



#' @rdname plot_umis_per_cell
#' @usage NULL
.plot_umis_per_cell_histogram <- function(object, min) {
    metrics <- metrics(object)
    ggplot(
        metrics,
        aes_(x = ~umi_counts,
             fill = ~sample_name)) +
        labs(title = "umi counts histogram",
             x = "umis per cell") +
        facet_wrap(~file_name) +
        geom_histogram(bins = bins) +
        geom_vline(alpha = 0.5,
                   color = warn_color,
                   size = 2L,
                   xintercept = min) +
        scale_x_log10() +
        scale_y_sqrt() +
        theme(axis.text.x = element_text(angle = 90L, hjust = 1L))
}



#' @rdname plot_umis_per_cell
#' @usage NULL
.plot_umis_per_cell <- function(object, min) {
    plot_grid(.plot_umis_per_cell_histogram(object, min),
              .plot_umis_per_cell_boxplot(object, min),
              labels = "auto",
              nrow = 2L)
}



#' @rdname plot_umis_per_cell
#' @export
setMethod(
    "plot_umis_per_cell",
    "bcbioSCDataSet",
    function(object, min = 1000L) {
        .plot_umis_per_cell(object, min)
    })

#' @rdname plot_umis_per_cell
#' @export
setMethod(
    "plot_umis_per_cell",
    "SCSubset",
    function(object) {
        min <- object %>%
            metadata %>%
            .[["filtering_criteria"]] %>%
            .[["min_umis"]]
        .plot_umis_per_cell(object, min)
    })
