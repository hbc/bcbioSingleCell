#' Counts Accessor
#'
#' @rdname counts
#' @name counts
#' @author Michael Steinbaugh
#'
#' @inheritParams all_generics
#' @param gene2symbol Convert Ensembl gene identifiers (rownames) to gene
#'   symbols. Recommended for passing counts to Seurat.
#' @param class Return counts as sparse matrix (**recommended**; `dgCMatrix`,
#'   `dgTMatrix`) or dense matrix (`matrix`).
#'
#' @return [matrix].
NULL



# Constructors ====
.counts <- function(
    object,
    gene2symbol = FALSE,
    class = "dgCMatrix") {
    supportedClasses <- c("dgCMatrix", "dgTMatrix", "matrix")
    if (!class %in% supportedClasses) {
        stop(paste("Supported classes:", toString(supportedClasses)))
    }
    counts <- assay(object)
    if (isTRUE(gene2symbol)) {
        rownames(counts) <- .gene2symbol(object)
    }
    as(counts, class)
}



# Methods ====
#' @rdname counts
#' @export
setMethod("counts", "bcbioSCDataSet", .counts)



#' @rdname counts
#' @export
setMethod("counts", "bcbioSCSubset", .counts)
