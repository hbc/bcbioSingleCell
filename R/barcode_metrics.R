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
    sparsecounts,
    metadata = NULL) {
    check_run(run)
    check_sparse(sparsecounts)

    # Ensembl annotations with broad class definitions
    annotations <- ensembl_annotations(run)

    coding <- annotations %>%
        filter_(.dots = quote(broad_class == "coding")) %>%
        .$external_gene_name %>% unique %>% sort
    mito <- annotations %>%
        filter_(.dots = quote(broad_class == "mito")) %>%
        .$external_gene_name %>% unique %>% sort

    # Calculate colSums
    metrics <- data.frame(
        identifier = colnames(sparsecounts),
        total_counts = Matrix::colSums(sparsecounts),
        genes_detected = Matrix::colSums(sparsecounts > 0),
        coding_counts = Matrix::colSums(sparsecounts[rownames(sparsecounts) %in% coding, ]),
        mito_counts = Matrix::colSums(sparsecounts[rownames(sparsecounts) %in% mito, ]))
    
    # dplyr operations
    metrics <- metrics %>%
        separate_("identifier",
                  c("sample_barcode", "cellular_barcode"),
                  sep = ":",
                  remove = FALSE) %>%
        mutate_(.dots = set_names(
            list(quote(log(genes_detected) / log(total_counts)),
                 quote(mito_counts / total_counts)),
            c("log_detected_per_count",
              "percent_mito")))

    # Attempt to detect metadata saved manually for consult, if necessary
    if (is.null(metadata)) {
        metafile <- list.files("meta",
                               pattern = "^sample_barcodes",
                               full.names = TRUE)
        if (length(metafile) == 1) {
            metadata <- read_sample_barcodes_metadata(metafile)
        } else if (length(metafile) > 1) {
            stop("multiple sample barcode metadata files matched...
                 need to clean the meta folder")
        } else {
            metadata <- NULL
        }
    }

    # Join metadata, if present
    if (!is.null(metadata)) {
        metrics <- left_join(metrics, metadata, by = "sample_barcode")
    }

    metrics <- metrics %>%
        arrange_(.dots = "identifier") %>%
        set_rownames(.$identifier)

    return(metrics)
}
