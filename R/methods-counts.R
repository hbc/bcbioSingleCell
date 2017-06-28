#' Counts accessor
#'
#' @rdname counts
#'
#' @author Michael Steinbaugh
#'
#' @param object [bcbioSCDataSet].
#' @param gene2symbol Convert Ensembl gene identifiers (rownames) to gene
#'   symbols. Recommended for passing counts to Seurat.
#' @param format Return counts as sparse matrix (**recommended**; `dgCMatrix`,
#'   `dgTMatrix`) or dense matrix (`matrix`).
#'
#' @return Counts matrix.



#' @rdname counts
#' @usage NULL
.counts <- function(
    object,
    gene2symbol = FALSE,
    format = "dgCMatrix") {
    supported_formats <- c("dgCMatrix", "dgTMatrix", "matrix")
    if (!format %in% supported_formats) {
        stop(paste("Supported formats", toString(supported_formats)))
    }
    counts <- assay(object)
    if (isTRUE(gene2symbol)) {
        message("Converting Ensembl gene identifiers to symbols")
        gene2symbol <- rowData(object) %>%
            as.data.frame %>%
            set_rownames(.[["ensgene"]]) %>%
            .[rownames(counts), c("ensgene", "symbol")]
        # Check for identifier degradation
        if (!identical(rownames(counts), gene2symbol[["ensgene"]]) |
            any(is.na(gene2symbol[["symbol"]]))) {
            stop(paste("Ensembl identifier degradation detected.",
                       "Please disable gene2symbol conversion."))
        }
        # Make gene symbols unique, if necessary
        symbol <- gene2symbol[["symbol"]]
        if (any(duplicated(symbol))) {
            symbol <- make.unique(symbol)
        }
        rownames(counts) <- symbol
    }
    as(counts, format)
}



#' @rdname counts
#' @export
setMethod("counts", "bcbioSCDataSet", .counts)

#' @rdname counts
#' @export
setMethod("counts", "SummarizedExperiment", .counts)
