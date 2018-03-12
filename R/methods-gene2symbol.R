#' Gene to Symbol Mappings
#'
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
#' # bcbioRNASeq ====
#' gene2symbol(bcb) %>% glimpse()
#'
#' # seurat ====
#' gene2symbol(seurat) %>% glimpse()
NULL



# Constructors =================================================================
.gene2symbolFromRowData <- function(object) {
    rowData <- rowData(object)
    if (is.null(rowData)) {
        abort("Object does not contain internal rowData")
    }
    cols <- c("geneID", "geneName")
    assert_is_subset(cols, colnames(rowData))
    rowData[, cols]
}



# Methods ======================================================================
#' @rdname gene2symbol
#' @export
setMethod(
    "gene2symbol",
    signature("seurat"),
    function(object) {
        gene2symbol <- bcbio(object, "gene2symbol")
        if (is.data.frame(gene2symbol)) {
            return(gene2symbol)
        }

        gene2symbol <- .gene2symbolFromRowData(object)
        if (is.data.frame(gene2symbol)) {
            return(gene2symbol)
        }

        rownames <- slot(object, "data") %>% rownames()
        if (!is.null(names(rownames))) {
            geneID <- names(rownames)
            geneName <- rownames
            gene2symbol <- data.frame(
                "geneID" = geneID,
                "geneName" = geneName,
                row.names = geneID,
                stringsAsFactors = FALSE
            )
            return(gene2symbol)
        }

        abort(paste(
            "seurat object doesn't appear to contain",
            "Ensembl gene-to-symbol mappings"
        ))
    }
)
