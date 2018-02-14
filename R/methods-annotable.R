#' Ensembl Annotations
#'
#' @rdname annotable
#' @name annotable
#' @author Michael Steinbaugh
#'
#' @importFrom basejump annotable
#'
#' @inheritParams AllGenerics
#'
#' @return [data.frame]
#'
#' @examples
#' load(system.file(
#'     file.path("extdata", "bcb.rda"),
#'     package = "bcbioSingleCell"))
#' load(system.file(
#'     file.path("extdata", "seurat.rda"),
#'     package = "bcbioSingleCell"))
#'
#' # bcbioSingleCell
#' annotable(bcb) %>% glimpse()
#'
#' # seurat
#' annotable(seurat) %>% glimpse()
NULL



# Constructors =================================================================
.annotable.bcbioSingleCell <- function(object) {  # nolint
    data <- rowData(object)
    rownames(data) <- slot(object, "NAMES")
    as.data.frame(data)
}



.annotable.seurat <- function(object) {  # nolint
    data <- bcbio(object, "annotable")
    assert_is_data.frame(data, severity = "warning")
    if (is.null(data)) {
        return(invisible())
    }
    # See if the gene identifiers are stashed as the names of the rownames
    genes <- slot(object, "data") %>%
        rownames() %>%
        names()
    assert_is_character(genes)
    data <- data[genes, , drop = FALSE]
    rownames(data) <- genes
    data
}



# Methods ======================================================================
#' @rdname annotable
#' @export
setMethod(
    "annotable",
    signature("bcbioSingleCell"),
    .annotable.bcbioSingleCell)



#' @rdname annotable
#' @export
setMethod(
    "annotable",
    signature("seurat"),
    .annotable.seurat)
