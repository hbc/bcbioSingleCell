#' Counts Accessor
#'
#' @rdname counts
#' @name counts
#' @author Michael Steinbaugh
#'
#' @inheritParams AllGenerics
#' @param gene2symbol Convert Ensembl gene identifiers (rownames) to gene
#'   symbols. Recommended for passing counts to Seurat.
#' @param genomeBuild *Optional*. If `gene2symbol = TRUE`, specify which
#'   genomeBuild to use for obtaining gene symbols.
#' @param as Return class (**recommended**; `dgCMatrix`,
#'   `dgTMatrix`) or dense matrix (`matrix`).
#'
#' @return Matrix class object.
NULL



# Constructors ====
.counts <- function(
    object,
    gene2symbol = FALSE,
    genomeBuild = NULL,
    as = "dgCMatrix") {
    supportedClasses <- c("dgCMatrix", "dgTMatrix", "matrix")
    if (!as %in% supportedClasses) {
        stop(paste("Supported classes:", toString(supportedClasses)),
             call. = FALSE)
    }
    counts <- assay(object)
    if (isTRUE(gene2symbol)) {
        counts <- gene2symbol(counts, genomeBuild = genomeBuild)
    }
    as(counts, as)
}



# Methods ====
#' @rdname counts
#' @export
setMethod("counts", "bcbioSCDataSet", .counts)

#' @rdname counts
#' @export
setMethod("counts", "bcbioSCSubset", .counts)
