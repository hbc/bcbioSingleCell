#' Generate barcode metrics summary
#'
#' @rdname metrics
#' @keywords internal
#' @author Rory Kirchner
#' @author Michael Steinbaugh
#'
#' @param run [bcbioSCDataSet].
#'
#' @return Tibble grouped by sample name.
.metrics <- function(run) {
    # Check for [Matrix::colSums()] methods support
    if (!"colSums,dgCMatrix-method" %in% methods(colSums)) {
        stop("dgCMatrix not supported in `colSums()`")
    }

    counts <- run$counts
    metadata <- run$metadata

    # Ensembl annotations, with broad class definitions (coding, mito)
    ensembl <- run$ensembl %>%
        ungroup %>%
        mutate(ensembl_transcript_id = NULL) %>%
        distinct
    coding <- ensembl %>%
        filter(.data$broad_class == "coding") %>%
        tidy_select(.data$external_gene_name) %>%
        .[[1]] %>% unique %>% sort
    mito <- ensembl %>%
        filter(.data$broad_class == "mito") %>%
        tidy_select(.data$external_gene_name) %>%
        .[[1]] %>% unique %>% sort

    data.frame(
        # Unique: `file_name` + `sample_barcode` + `cellular_barcode`
        unique = colnames(counts),
        total_counts = Matrix::colSums(counts),
        genes_detected = Matrix::colSums(counts > 0),
        coding_counts = Matrix::colSums(
            counts[rownames(counts) %in% coding, ]),
        mito_counts = Matrix::colSums(
            counts[rownames(counts) %in% mito, ])) %>%
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
        # Select sample name first
        tidy_select(.data$sample_name,
                    .data$sample_barcode,
                    .data$cellular_barcode,
                    everything()) %>%
        group_by(!!!syms(c("sample_name",  "sample_barcode"))) %>%
        arrange(desc(!!sym("total_counts")), .by_group = TRUE) %>%
        ungroup %>%
        mutate(rowname = paste(.data$sample_barcode,
                               .data$cellular_barcode,
                               sep = ":")) %>%
        as.data.frame %>%
        column_to_rownames %>%
        DataFrame
}











#' Unlist cellular barcodes
#'
#' Convert named list of cellular barcodes per sample to a data frame.
#'
#' @rdname unlist_barcodes
#' @keywords internal
#'
#' @author Michael Steinbaugh
#'
#' @param run [bcbioSCDataSet].
.unlist_barcodes <- function(run) {
    barcodes <- run$barcodes
    message("Converting nested barcodes to data frame...")
    pbmclapply(seq_along(barcodes), function(a) {
        barcodes[a] %>%
            as.data.frame %>%
            rownames_to_column %>%
            set_names(c("cellular_barcode", "reads")) %>%
            arrange(!!!syms(c("reads", "cellular_barcode"))) %>%
            mutate(log10_reads = log10(.data$reads),
                   sample_barcode = names(barcodes[a]))
    }) %>% bind_rows %>%
        left_join(run$metadata[, c("sample_barcode", "sample_name")],
                  by = "sample_barcode")
}
