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
#'
#' @return Data frame
#' @export
barcode_metrics <- function(counts,
                            annotations,
                            metadata=NULL) {
    coding <- annotations %>%
        dplyr::filter_(.dots = quote(broad_class == "coding")) %>%
        .$gene_name %>% unique %>% sort
    mito <- annotations %>%
        dplyr::filter_(.dots = quote(broad_class == "mito")) %>%
        .$gene_name %>% unique %>% sort

    # `rmarkdown::render()` doesn't handle `dgTMatrix` objects properly.
    # `colSums(counts)` fails here unless we coerce `counts` to a matrix first.
    counts <- as.matrix(counts)

    metrics <- data.frame(
        identifier = colnames(counts),
        total_counts = colSums(counts),
        genes_detected = colSums(counts > 0),
        coding_counts = colSums(counts[rownames(counts) %in% coding, ]),
        mito_counts = colSums(counts[rownames(counts) %in% mito, ])
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
            ))
    if(!is.null(metadata)) {
      metrics <- metrics %>%
        dplyr::left_join(metadata[, c("sample_barcode", "sample")],
                         by = "sample_barcode")
    }
    metrics <- metrics %>%
      dplyr::arrange_(.dots = "identifier") %>%
      set_rownames("identifier")

    return(metrics)
}
