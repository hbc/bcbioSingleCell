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



# Methods ====
#' @rdname annotable
#' @export
setMethod(
    "annotable",
    signature("bcbioSingleCell"),
    function(object) {
        annotable <- as.data.frame(rowData(object))
        rownames(annotable) <- slot(object, "NAMES")
        annotable
    })



#' @rdname annotable
#' @export
setMethod(
    "annotable",
    signature("seurat"),
    function(object) {
        annotable <- bcbio(object, "annotable")
        if (is.null(annotable)) return(NULL)
        rownames <- slot(object, "data") %>% rownames() %>% names()
        if (is.null(rownames)) return(NULL)
        annotable <- annotable[rownames, ]
        rownames(annotable) <- rownames
        annotable
    })
