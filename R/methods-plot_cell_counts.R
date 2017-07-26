#' Plot Cell Counts
#'
#' @rdname plot_cell_counts
#' @author Michael Steinbaugh, Rory Kirchner
#' @inherit plot_genes_per_cell



#' @rdname plot_cell_counts
#' @usage NULL
.plot_cell_counts <- function(object) {
    metrics <- metrics(object)
    cell_counts <- metrics %>%
        group_by(!!sym("sample_id")) %>%
        summarise(cells = n()) %>%
        left_join(sample_metadata(object), by = "sample_id")
    interesting_group <- interesting_groups(object)[[1L]]
    p <- ggplot(cell_counts,
                aes_(x = ~sample_name,
                     y = ~cells,
                     fill = as.name(interesting_group))) +
        labs(x = "sample",
             y = "cell count") +
        geom_bar(stat = "identity") +
        geom_text(vjust = -0.5, aes_(label = ~cells)) +
        theme(axis.text.x = element_text(angle = 90L, hjust = 1L))
    if (isTRUE(metadata(object)[["multiplexed_fastq"]])) {
        p <- p + facet_wrap(~file_name)
    }
    p
}



#' @rdname plot_cell_counts
#' @export
setMethod("plot_cell_counts", "bcbioSCDataSet", .plot_cell_counts)



#' @rdname plot_cell_counts
#' @export
setMethod("plot_cell_counts", "bcbioSCSubset", .plot_cell_counts)
