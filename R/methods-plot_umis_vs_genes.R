#' Plot UMI and Gene Correlation
#'
#' @rdname plot_umis_vs_genes
#' @author Michael Steinbaugh, Rory Kirchner
#' @inherit plot_genes_per_cell



#' @rdname plot_umis_vs_genes
#' @usage NULL
.plot_umis_vs_genes <- function(object) {
    metrics <- metrics(object)
    p <- ggplot(metrics,
        aes_(x = ~umi_counts,
             y = ~genes_detected,
             color = ~sample_name)) +
        labs(x = "umis per cell",
             y = "genes per cell") +
        geom_smooth(method = "lm", se = FALSE) +
        scale_x_log10() +
        scale_y_log10() +
        theme(axis.text.x = element_text(angle = 90L, hjust = 1L))
    if (isTRUE(metadata(object)[["multiplexed_fastq"]])) {
        p <- p + facet_wrap(~file_name)
    }
    p
}



#' @rdname plot_umis_vs_genes
#' @export
setMethod("plot_umis_vs_genes", "bcbioSCDataSet", .plot_umis_vs_genes)



#' @rdname plot_umis_vs_genes
#' @export
setMethod("plot_umis_vs_genes", "bcbioSCSubset", .plot_umis_vs_genes)
