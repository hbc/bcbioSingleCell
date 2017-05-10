#' Generate barcode metrics summary
#'
#' @author Rory Kirchner
#' @author Michael Steinbaugh
#'
#' @param run \code{bcbio-nextgen} run
#' @param counts Sparse counts matrix
#'
#' @return Data frame
#' @export
barcode_metrics <- function(
    run,
    counts) {
    check_run(run)
    check_sparse(counts)
    ensembl <- run$ensembl
    metadata <- run$metadata

    coding <- ensembl %>%
        .[.$broad_class == "coding", ] %>%
        .$external_gene_name %>% unique %>% sort
    mito <- ensembl %>%
        .[.$broad_class == "mito", ] %>%
        .$external_gene_name %>% unique %>% sort

    metrics <- tibble(
        # Unique: `file_name` + `sample_barcode` + `cellular_barcode`
        unique = colnames(counts),
        # Matrix::colSums()
        total_counts = colSums(counts),
        genes_detected = colSums(counts > 0),
        coding_counts = colSums(
            counts[rownames(counts) %in% coding, ]),
        mito_counts = colSums(
            counts[rownames(counts) %in% mito, ])
    ) %>%
        # Separate the barcodes, later used to join metadata
        separate_("unique",
                  c("sample_barcode", "cellular_barcode"),
                  sep = ":",
                  remove = FALSE) %>%
        # Summary statistics
        mutate(log10_detected_per_count = log10(.data$genes_detected) /
                   log10(.data$total_counts),
               percent_mito = .data$mito_counts / .data$total_counts) %>%
        # Join the sample names
        left_join(metadata[, c("sample_barcode", "sample_name")],
                  by = "sample_barcode") %>%
        # Arrange by unique identifier
        # arrange(!!sym("unique")) %>%
        .[order(.$unique), ] %>%
        # Filter barcodes that don't match samples of interest
        # filter(!is.na(.data$sample_name)) %>%
        .[!is.na(.$sample_name), ] %>%
        as.data.frame %>%
        set_rownames(.$unique)

    return(metrics)
}
