#' Aggregate Features
#'
#' @rdname aggregateFeatures
#' @name aggregateFeatures
#' @family Data Management Utilities
#' @author Rory Kirchner, Michael Steinbaugh
#' @keywords internal
#'
#' @inheritParams AllGenerics
#'
#' @param object Sparse counts matrix (e.g. [dgCMatrix]).
#' @param featureids Feature identifiers (e.g. gene or transcript IDs).
#'
#' @return [dgCMatrix].
#'
#' @return [bcbioSingleCell].
NULL



# Constructors ====
#' @importFrom Matrix.utils aggregate.Matrix
.aggregateSparseFeatures <- function(object, featureids) {
    rownames(sparse) <- featureids
    sparse <- sparse[!is.na(rownames(sparse)), , drop = FALSE]
    aggregate.Matrix(sparse, groupings = rownames(sparse), fun = "sum")
}



# Methods ====
#' @rdname aggregateFeatures
#' @export
setMethod(
    "aggregateFeatures",
    signature("bcbioSingleCell"),
    function(object, featureids) {
        warning(paste(
            "Draft function.",
            "Returning an aggregated counts matrix."),
            call. = FALSE)
        .aggregateSparseFeatures(
            object = assay(object),
            featureids = featureids
        )
    })



#' @rdname aggregateFeatures
#' @export
setMethod(
    "aggregateFeatures",
    signature("dgCMatrix"),
    .aggregateSparseFeatures)



#' @rdname aggregateFeatures
#' @export
setMethod(
    "aggregateFeatures",
    signature("dgTMatrix"),
    .aggregateSparseFeatures)
