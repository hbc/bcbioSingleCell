#' Filter cellular barcodes
#'
#' Apply gene detection, mitochondrial abundance, and novelty score cutoffs to
#' cellular barcodes.
#'
#' @author Michael Steinbaugh
#'
#' @param metrics Barcode metrics data frame
#' @param min_genes Minimum number of genes detected
#' @param max_genes Maximum number of genes detected
#' @param percent_mito Maximum mitochondrial pct abundance (\code{0-1} scale)
#' @param novelty Minimum novelty score
#' @param plot Print relevant plots
#'
#' @return Filtered metrics data frame
#' @export
filter_barcodes <- function(
    metrics,
    min_genes = 500,
    max_genes = 5000,
    percent_mito = 0.25,
    novelty = 0.75,
    plot = TRUE) {
    filtered <- metrics %>%
        .[.$genes_detected > min_genes, ] %>%
        .[.$genes_detected < max_genes, ] %>%
        .[.$percent_mito < percent_mito, ] %>%
        .[.$log10_detected_per_count > novelty, ]
    if (isTRUE(plot)) {
        show(plot_total_cells(filtered))
        plot_total_counts(filtered)
        plot_genes_detected(filtered)
        show(plot_total_vs_detected(filtered))
        plot_mito_counts(filtered)
        plot_novelty(filtered)
    }
    return(filtered)
}
