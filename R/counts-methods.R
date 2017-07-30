#' Counts Accessor
#'
#' @rdname counts
#' @name counts
#' @author Michael Steinbaugh
#'
#' @inheritParams AllGenerics
#'
#' @param gene2symbol Convert Ensembl gene identifiers (rownames) to gene
#'   symbols. Recommended for passing counts to Seurat.
#' @param class Return counts as sparse matrix (**recommended**; `dgCMatrix`,
#'   `dgTMatrix`) or dense matrix (`matrix`).
#'
#' @return Matrix class object.
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
        g2s <- object %>%
            rownames %>%
            .[[1L]] %>%
            detectOrganism %>%
            gene2symbol %>%
            .[rownames(counts), ]
        rownames(counts) <- g2s[rownames(counts), "symbol"] %>% make.unique
        # Warn if any symbols are NA
        if (any(is.na(rownames(counts)))) {
            warning("NA symbols detected in counts matrix")
        }
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
