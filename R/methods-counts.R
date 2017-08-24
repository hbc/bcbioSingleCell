#' Counts Accessor
#'
#' @rdname counts
#' @name counts
#' @author Michael Steinbaugh
#'
#' @inheritParams AllGenerics
#' @param gene2symbol Convert Ensembl gene identifiers (rownames) to gene
#'   symbols. Recommended for passing counts to Seurat.
#' @param as Return class (**recommended**; `dgCMatrix`,
#'   `dgTMatrix`) or dense matrix (`matrix`).
#'
#' @return Matrix class object.
NULL



# Constructors ====
.counts <- function(
    object,
    gene2symbol = FALSE,
    as = "dgCMatrix") {
    supportedClasses <- c("dgCMatrix", "dgTMatrix", "matrix")
    if (!as %in% supportedClasses) {
        stop(paste("Supported classes:", toString(supportedClasses)),
             call. = FALSE)
    }
    counts <- assay(object)
    if (isTRUE(gene2symbol)) {
        g2s <- metadata(object)[["gene2symbol"]]
        if (!all(rownames(counts) %in% rownames(g2s))) {
            rescale <- g2s %>%
                .[rownames(counts), ] %>%
                set_rownames(rownames(counts))
            matched <- rescale %>%
                .[!is.na(.[["symbol"]]), ]
            unmatched <- rescale %>%
                .[is.na(.[["symbol"]]), ]
            warning(paste(
                "Unmatched in gene2symbol:",
                toString(rownames(unmatched))), call. = FALSE)
            unmatched[["ensgene"]] <- rownames(unmatched)
            unmatched[["symbol"]] <- rownames(unmatched)
            g2s <- rbind(matched, unmatched)
        }
        g2s <- g2s[rownames(counts), ]
        rows <- pull(g2s, "symbol")
        names(rows) <- rownames(g2s)
        rownames(counts) <- rows
    }
    as(counts, as)
}



# Methods ====
#' @rdname counts
#' @export
setMethod("counts", "bcbioSCDataSet", .counts)



#' @rdname counts
#' @export
setMethod("counts", "bcbioSCFiltered", .counts)
