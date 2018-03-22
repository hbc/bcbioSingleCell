# FIXME Rename `gene2symbol` to `convertGenesToSymbols. Catch the legacy param.

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
#' @param convertGenesToSymbols *Not recommended.* Convert Ensembl gene identifiers
#'   (rownames) to gene symbols. Required for passing counts to Seurat.
#' @param normalized Normalized (`TRUE`) or raw (`FALSE`) counts.
#'
#' @return Counts matrix.
#'
#' @examples
#' load(system.file("extdata/bcb.rda", package = "bcbioSingleCell"))
#' load(system.file("extdata/seurat.rda", package = "bcbioSingleCell"))
#'
#' # bcbioSingleCell
#' counts(bcb, convertGenesToSymbols = FALSE) %>% glimpse()
#' counts(bcb, convertGenesToSymbols = TRUE) %>% glimpse()
#'
#' # seurat
#' counts(seurat, normalized = FALSE) %>% glimpse()
#' counts(seurat, normalized = TRUE) %>% glimpse()
NULL



# Constructors =================================================================
.counts <- function(
    object,
    convertGenesToSymbols = FALSE
) {
    counts <- assay(object)
    if (isTRUE(convertGenesToSymbols)) {
        g2s <- gene2symbol(object)
        assert_is_subset(rownames(counts), rownames(g2s))
        # Resize the gene2symbol data.frame
        g2s <- g2s[rownames(counts), , drop = FALSE]
        # Abort on any NA gene names
        stopifnot(!any(is.na(g2s[["geneName"]])))
        rownames <- make.names(g2s[["geneName"]], unique = TRUE)
        names(rownames) <- g2s[["geneID"]]
        rownames(counts) <- rownames
    }
    counts
}



# Methods ======================================================================
#' @rdname counts
#' @export
setMethod(
    "counts",
    signature("bcbioSingleCell"),
    .counts
)



#' @rdname counts
#' @export
setMethod(
    "counts",
    signature("seurat"),
    function(object, normalized = FALSE) {
        assert_is_a_bool(normalized)
        # seurat also stashes scaled counts in `scale.data`
        if (identical(normalized, FALSE)) {
            slot(object, "raw.data")
        } else if (identical(normalized, TRUE)) {
            slot(object, "data")
        }
    }
)
