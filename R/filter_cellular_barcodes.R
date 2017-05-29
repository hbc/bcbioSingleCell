#' Filter cellular barcodes
#'
#' Apply gene detection, mitochondrial abundance, and novelty score cutoffs to
#' cellular barcodes.
#'
#' @author Michael Steinbaugh
#'
#' @param run [bcbioSCDataSet].
#' @param min_genes Minimum number of genes detected.
#' @param max_genes Maximum number of genes detected.
#' @param mito_ratio Maximum relative mitochondrial abundance (`0-1` scale).
#' @param novelty Minimum novelty score.
#' @param plot Print relevant plots.
#'
#' @return Filtered [bcbioSCDataSet].
#' @export
filter_cellular_barcodes <- function(
    run,
    min_genes = get("min_genes", envir = parent.frame()),
    max_genes = get("max_genes", envir = parent.frame()),
    mito_ratio = get("mito_ratio", envir = parent.frame()),
    novelty = get("novelty", envir = parent.frame()),
    plot = TRUE) {
    check_run(run)
    run$metrics <- run$metrics %>%
        filter(.data$genes_detected > !!min_genes,
               .data$genes_detected < !!max_genes,
               .data$mito_ratio < !!mito_ratio,
               .data$log10_detected_per_count > !!novelty)
    run$filtered <- TRUE
    if (isTRUE(plot)) {
        show(plot_total_cells(run))
        plot_total_counts(run)
        plot_genes_detected(run, min_genes = min_genes, max_genes = max_genes)
        show(plot_total_vs_detected(run))
        plot_mito(run, mito_ratio = mito_ratio)
        plot_novelty(run, novelty = novelty)
    }
    run
}



# Globals used for QC plot labels ====
#' @rdname filter_cellular_barcodes
#' @description Default minimum gene count cutoff.
#' @export
#' @examples
#' min_genes
min_genes <- 500

#' @rdname filter_cellular_barcodes
#' @description Default maximum gene count cutoff.
#' @export
#' @examples
#' max_genes
max_genes <- 5000

#' @rdname filter_cellular_barcodes
#' @description Default maximum relative mitochondrial abundance.
#' @export
#' @examples
#' mito_ratio
mito_ratio <- 0.25

#' @rdname filter_cellular_barcodes
#' @description Default minimum novelty score.
#' @export
#' @examples
#' novelty
novelty <- 0.75
