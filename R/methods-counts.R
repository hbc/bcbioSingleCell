#' Counts Accessor
#'
#' @rdname counts
#' @inheritParams all_generics
#'
#' @author Michael Steinbaugh
#'
#' @param gene2symbol Convert Ensembl gene identifiers (rownames) to gene
#'   symbols. Recommended for passing counts to Seurat.
#' @param format Return counts as sparse matrix (**recommended**; `dgCMatrix`,
#'   `dgTMatrix`) or dense matrix (`matrix`).
#'
#' @return [matrix].



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
        rownames(counts) <- .gene2symbol(object)
    }
    as(counts, format)
}



#' @rdname counts
#' @export
setMethod("counts", "bcbioSCDataSet", .counts)

#' @rdname counts
#' @export
setMethod("counts", "bcbioSCSubset", .counts)
