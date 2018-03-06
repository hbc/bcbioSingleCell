#' Gene to Symbol Mappings
#'
#' @rdname gene2symbol
#' @name gene2symbol
#' @author Michael Steinbaugh
#'
#' @importFrom basejump gene2symbol
#'
#' @inheritParams general
#'
#' @examples
#' load(system.file("extdata/bcb.rda", package = "bcbioSingleCell"))
#' load(system.file("extdata/seurat.rda", package = "bcbioSingleCell"))
#'
#' # bcbioRNASeq
#' gene2symbol(bcb) %>% glimpse()
#'
#' # seurat
#' gene2symbol(seurat) %>% glimpse()
NULL



# Constructors =================================================================
.gene2symbolFromRowData <- function(object) {
    rowData <- rowData(object)
    if (is.null(rowData)) {
        abort("Object does not contain internal rowData")
    }
    cols <- c("ensgene", "symbol")
    if (!all(cols %in% colnames(rowData))) {
        abort(paste(
            toString(cols),
            "missing from internal rowData"
        ))
    }
    rowData[, cols]
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

        # Fall back to generating from rowData
        .gene2symbolFromRowData(object)
    })



#' @rdname gene2symbol
#' @export
setMethod(
    "gene2symbol",
    signature("seurat"),
    function(object) {
        gene2symbol <- bcbio(object, "gene2symbol")
        if (is.data.frame(gene2symbol)) return(gene2symbol)

        gene2symbol <- .gene2symbolFromRowData(object)
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

        warn(paste(
            "seurat object doesn't appear to contain",
            "Ensembl gene to symbol mappings"
        ))
    }
)
