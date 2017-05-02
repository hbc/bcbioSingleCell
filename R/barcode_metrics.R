#' Generate barcode metrics summary
#'
#' @author Rory Kirchner
#' @author Michael Steinbaugh
#'
#' @param run \code{bcbio-nextgen} run
#' @param sparsecounts Sparse counts matrix
#' @param metadata Sample barcodes metadata
#'
#' @return Data frame
#' @export
barcode_metrics <- function(
    run,
    sparsecounts) {
    check_run(run)
    check_sparse(sparsecounts)
    ensembl <- run$ensembl
    metadata <- run$metadata

    coding <- ensembl %>%
        .[.$broad_class == "coding", ] %>%
        .$external_gene_name %>% unique %>% sort
    mito <- ensembl %>%
        .[.$broad_class == "mito", ] %>%
        .$external_gene_name %>% unique %>% sort

    # Calculate colSums
    metrics <- tibble(
        unique = colnames(sparsecounts),
        total_counts = Matrix::colSums(sparsecounts),
        genes_detected = Matrix::colSums(sparsecounts > 0),
        coding_counts = Matrix::colSums(
            sparsecounts[rownames(sparsecounts) %in% coding, ]),
        mito_counts = Matrix::colSums(
            sparsecounts[rownames(sparsecounts) %in% mito, ])) %>%
        # Unique identifier
        separate_("unique",
                  c("sample_barcode", "cellular_barcode"),
                  sep = ":",
                  remove = FALSE) %>%
        # Summary statistics
        mutate(log10_detected_per_count = log10(.data$genes_detected) /
                   log10(.data$total_counts),
               percent_mito = .data$mito_counts / .data$total_counts)

    # Join metadata, if present
    if (!is.null(metadata)) {
        metrics <- left_join(metrics, metadata, by = "sample_barcode")
    }

    metrics <- metrics[order(metrics$unique), ]
    return(metrics)
}
