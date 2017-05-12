#' Filter cellular barcodes.
#'
#' Apply gene detection, mitochondrial abundance, and novelty score cutoffs to
#' cellular barcodes.
#'
#' @author Michael Steinbaugh
#'
#' @param run bcbio-nextgen scRNA-seq run.
#' @param min_genes Minimum number of genes detected.
#' @param max_genes Maximum number of genes detected.
#' @param percent_mito Maximum mitochondrial pct abundance (\code{0-1} scale).
#' @param novelty Minimum novelty score.
#' @param plot Print relevant plots.
#'
#' @return Filtered metrics data frame.
#' @export
filter_cellular_barcodes <- function(
    run,
    min_genes = 500,
    max_genes = 5000,
    percent_mito = 25,
    novelty = 0.75,
    plot = TRUE) {
    check_run(run)
    filtered <- run$metrics %>%
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



# Globals used for QC plot labels ====
#' @rdname filter_cellular_barcodes
#' @description Default minimum gene count cutoff (500)
#' @export
#' @examples
#' min_genes
min_genes <- 500

#' @rdname filter_cellular_barcodes
#' @description Default maximum gene count cutoff (5000)
#' @export
#' @examples
#' max_genes
max_genes <- 5000

#' @rdname filter_cellular_barcodes
#' @description Default minimum novelty score (0.75)
#' @export
#' @examples
#' novelty
novelty <- 0.75

#' @rdname filter_cellular_barcodes
#' @description Default maximum mitochondrial percent abundance (25)
#' @export
#' @examples
#' percent_mito
percent_mito <- 25
