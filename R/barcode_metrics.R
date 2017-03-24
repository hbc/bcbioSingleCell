#' Generate barcode metrics summary
#'
#' @author Rory Kirchner
#' @author Michael Steinbaugh
#'
#' @import dplyr
#' @import tidyr
#'
#' @param counts `bcbio-nextgen` scRNA-seq counts sparse matrix
#' @param annotations Ensembl annotations data frame
#' @param metadata Sample metadata data frame
#' @param tx2gene Convert transcript-level annotations to gene-level
#'
#' @return Data frame
#' @export
barcode_metrics <- function(counts,
                            annotations,
                            metadata,
                            tx2gene = TRUE) {
    if (isTRUE(tx2gene)) {
        annotations$ensembl_transcript_id <- NULL
        annotations <- annotations %>%
            dplyr::arrange_(.dots = "ensembl_gene_id") %>%
            dplyr::distinct(.) %>%
            set_rownames("ensembl_gene_id")
    }

    coding <- annotations %>%
        dplyr::filter_(.dots = quote(broad_class == "coding")) %>%
        .$ensembl_transcript_id %>% unique %>% sort
    mito <- annotations %>%
        dplyr::filter_(.dots = quote(broad_class == "mito")) %>%
        .$ensembl_transcript_id %>% unique %>% sort

    # `rmarkdown::render()` doesn't handle `dgTMatrix` objects properly.
    # `colSums(counts)` fails here unless we coerce `counts` to a matrix first.
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
        dplyr::mutate_(.dots = set_names(
            list(quote(log(genes_detected) / log(total_counts)),
                 quote(mito_counts / total_counts)),
            c("log_detected_per_count",
              "percent_mito")
        )) %>%
        dplyr::arrange_(.dots = "identifier") %>%
        dplyr::left_join(metadata[, c("sample_barcode", "sample")],
                         by = "sample_barcode") %>%
        set_rownames("identifier")

    return(metrics)
}
