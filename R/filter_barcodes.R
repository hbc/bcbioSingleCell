#' Filter cellular barcodes
#'
#' Apply gene detection, mitochondrial abundance, and novelty score cutoffs to
#' cellular barcodes.
#'
#' @author Michael Steinbaugh
#'
#' @import dplyr
#' @importFrom basejump setRownames
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
filter_barcodes <- function(metrics,
                            min_genes = 500,
                            max_genes = 5000,
                            percent_mito = 0.2,
                            novelty = 0.8,
                            plot = TRUE) {
    filtered <- metrics %>%
        dplyr::filter_(.dots = list(
            ~genes_detected > min_genes,
            ~genes_detected < max_genes,
            # Names need to be distinct, right?
            ~percent_mito < get("percent_mito", envir = environment()),
            ~log_detected_per_count > novelty
        )) %>%
        basejump::setRownames(., "identifier")

    if (isTRUE(plot)) {
        plot_total_cells(filtered)
        plot_genes_detected(filtered)
        plot_total_vs_detected(filtered)
    }

    return(filtered)
}
