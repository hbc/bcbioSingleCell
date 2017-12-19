#' Gene to Symbol Mappings
#'
#' @rdname gene2symbol
#' @name gene2symbol
#' @author Michael Steinbaugh
#'
#' @importFrom basejump gene2symbol
#'
#' @inheritParams AllGenerics
#'
#' @examples
#' load(system.file(
#'     file.path("extdata", "bcb.rda"),
#'     package = "bcbioSingleCell"))
#' load(system.file(
#'     file.path("extdata", "seurat.rda"),
#'     package = "bcbioSingleCell"))
#'
#' # bcbioRNASeq
#' gene2symbol(bcb) %>% glimpse()
#'
#' # seurat
#' gene2symbol(seurat) %>% glimpse()
NULL



# Constructors =================================================================
.gene2symbolFromAnnotable <- function(object) {
    annotable <- annotable(object)
    if (is.null(annotable)) {
        stop("Object does not contain internal annotable")
    }
    cols <- c("ensgene", "symbol")
    if (!all(cols %in% colnames(annotable))) {
        stop(paste(
            toString(cols),
            "missing from internal annotable"
        ), call. = FALSE)
    }
    annotable[, cols]
}


# Methods ======================================================================
#' @rdname gene2symbol
#' @export
setMethod(
    "gene2symbol",
    signature("bcbioSingleCell"),
    function(object) {
        gene2symbol <- metadata(object)[["gene2symbol"]]
        if (is.data.frame(gene2symbol)) return(gene2symbol)

        # Fall back to generating from annotable
        .gene2symbolFromAnnotable(object)
    })



#' @rdname gene2symbol
#' @export
setMethod(
    "gene2symbol",
    signature("seurat"),
    function(object) {
        gene2symbol <- bcbio(object, "gene2symbol")
        if (is.data.frame(gene2symbol)) return(gene2symbol)

        gene2symbol <- .gene2symbolFromAnnotable(object)
        if (is.data.frame(gene2symbol)) return(gene2symbol)

        rownames <- slot(object, "data") %>% rownames()
        if (!is.null(names(rownames))) {
            ensgene <- names(rownames)
            symbol <- rownames
            gene2symbol <- data.frame(
                "ensgene" = ensgene,
                "symbol" = symbol,
                row.names = ensgene,
                stringsAsFactors = FALSE
            )
            return(gene2symbol)
        }

        warning(paste(
            "seurat object doesn't appear to contain",
            "Ensembl gene to symbol mappings"
        ))
    })
