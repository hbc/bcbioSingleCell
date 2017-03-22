#' Generate barcode metrics summary
#'
#' @author Rory Kirchner
#' @author Michael Steinbaugh
#'
#' @import dplyr
#' @import tidyr
#' @importFrom basejump setRownames
#' @importFrom stats setNames
#'
#' @param counts `bcbio-nextgen` scRNA-seq counts sparse matrix
#' @param annotations Ensembl annotations data frame
#' @param metadata Sample metadata data frame
#'
#' @return Data frame
#' @export
barcode_metrics <- function(counts, annotations, metadata) {
    coding <- annotations %>%
        # dplyr::filter(broad_class == "coding") %>%
        dplyr::filter_(.dots = quote(broad_class == "coding")) %>%
        .$ensembl_transcript_id %>% unique %>% sort
    mito <- annotations %>%
        dplyr::filter_(.dots = quote(broad_class == "mito")) %>%
        .$ensembl_transcript_id %>% unique %>% sort

    # `rmarkdown::render()` handle `dgTMatrix` objects properly.
    # Fails on `colSums(counts)` unless we coerce `counts` to a matrix first.
    counts_matrix <- as.matrix(counts)

    metrics <- data.frame(
        identifier = colnames(counts_matrix),
        total_counts = colSums(counts_matrix),
        genes_detected = colSums(counts_matrix > 0),
        coding_counts = colSums(counts_matrix[rownames(counts_matrix) %in% coding, ]),
        mito_counts = colSums(counts_matrix[rownames(counts_matrix) %in% mito, ])
    ) %>%
        tidyr::separate_("identifier",
                         c("sample_barcode", "cellular_barcode"),
                         sep = ":",
                         remove = FALSE) %>%
        dplyr::mutate_(.dots = stats::setNames(
            list(quote(log(genes_detected) / log(total_counts)),
                 quote(mito_counts / total_counts)),
            c("log_detected_per_count",
              "percent_mito")
        )) %>%
        dplyr::arrange_(.dots = "identifier") %>%
        dplyr::left_join(metadata[, c("sample_barcode", "sample")],
                         by = "sample_barcode") %>%
        basejump::setRownames(., "identifier")

    return(metrics)
}
