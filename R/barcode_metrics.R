#' Generate barcode metrics summary.
#'
#' @author Rory Kirchner
#' @author Michael Steinbaugh
#'
#' @param run bcbio-nextgen run.
#'
#' @return Tibble grouped by sample name.
#' @export
barcode_metrics <- function(run) {
    import_tidy_verbs()
    check_run(run)
    counts <- run$counts
    ensembl <- run$ensembl
    metadata <- run$metadata

    coding <- ensembl %>%
        filter(.data$broad_class == "coding") %>%
        select(.data$external_gene_name) %>%
        .[[1]] %>% unique %>% sort
    mito <- ensembl %>%
        filter(.data$broad_class == "mito") %>%
        select(.data$external_gene_name) %>%
        .[[1]] %>% unique %>% sort

    metrics <- tibble(
        # Unique: `file_name` + `sample_barcode` + `cellular_barcode`
        unique = colnames(counts),
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
        mutate(unique = NULL,
               log10_detected_per_count = log10(.data$genes_detected) /
                   log10(.data$total_counts),
               mito_ratio = .data$mito_counts / .data$total_counts) %>%
        left_join(metadata[, c("sample_barcode", "sample_name")],
                  by = "sample_barcode") %>%
        # Place sample name first
        select(.data$sample_name, everything()) %>%
        # Filter barcodes matching samples
        filter(!is.na(.data$sample_name)) %>%
        group_by(!!sym("sample_name")) %>%
        arrange(!!sym("sample_name"), .by_group = TRUE)

    return(metrics)
}
