#' Counts Accessor
#'
#' @rdname counts
#' @name counts
#' @author Michael Steinbaugh
#'
#' @importFrom BiocGenerics counts
#'
#' @inheritParams general
#'
#' @param gene2symbol Convert Ensembl gene identifiers (rownames) to gene
#'   symbols. Recommended for passing counts to Seurat.
#' @param normalized Normalized (`TRUE`) or raw (`FALSE`) counts.
#'
#' @return Counts matrix.
#'
#' @examples
#' load(system.file("extdata/bcb.rda", package = "bcbioSingleCell"))
#' load(system.file("extdata/seurat.rda", package = "bcbioSingleCell"))
#'
#' # bcbioSingleCell
#' counts(bcb, gene2symbol = FALSE) %>% glimpse()
#' counts(bcb, gene2symbol = TRUE) %>% glimpse()
#'
#' # seurat
#' counts(seurat, normalized = FALSE) %>% glimpse()
#' counts(seurat, normalized = TRUE) %>% glimpse()
NULL



# Constructors =================================================================
#' @importFrom dplyr pull
#' @importFrom magrittr set_rownames
.counts <- function(
    object,
    gene2symbol = FALSE) {
    counts <- assay(object)
    if (isTRUE(gene2symbol)) {
        g2s <- metadata(object)[["gene2symbol"]]
        if (!all(rownames(counts) %in% rownames(g2s))) {
            # Resize the gene2symbol data.frame
            g2s <- g2s %>%
                .[rownames(counts), , drop = FALSE] %>%
                set_rownames(rownames(counts))
            matched <- g2s %>%
                .[!is.na(.[["geneName"]]), , drop = FALSE]
            unmatched <- g2s %>%
                .[is.na(.[["geneName"]]), , drop = FALSE]
            warn(paste(
                "Unmatched in gene2symbol:",
                toString(rownames(unmatched))
            ))
            unmatched[["geneID"]] <- rownames(unmatched)
            unmatched[["geneName"]] <- rownames(unmatched)
            g2s <- rbind(matched, unmatched)
        }
        g2s <- g2s[rownames(counts), , drop = FALSE]
        rows <- pull(g2s, "geneName")
        names(rows) <- rownames(g2s)
        rownames(counts) <- rows
    }
    counts
}



# Methods ======================================================================
#' @rdname counts
#' @export
setMethod(
    "counts",
    signature("bcbioSingleCell"),
    .counts)



#' @rdname counts
#' @export
setMethod(
    "counts",
    signature("seurat"),
    function(object, normalized = FALSE) {
        if (identical(normalized, FALSE)) {
            slot(object, "raw.data")
        } else if (identical(normalized, TRUE)) {
            slot(object, "data")
        } else if (normalized == "scaled") {
            slot(object, "scale.data")
        } else {
            warn("Unsupported `normalized` argument")
            return(NULL)
        }
    }
)
