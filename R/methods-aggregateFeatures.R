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
#' @return [bcbioSingleCell].
NULL



# Constructors ====
#' Aggregate Features Constructor
#'
#' @author Rory Kirchner
#'
#' @importFrom Matrix.utils aggregate.Matrix
#'
#' @param sparse Sparse counts matrix (e.g. [dgCMatrix]).
#' @param featureids Feature identifiers (e.g. gene or transcript IDs).
#'
#' @return [dgCMatrix].
#' @noRd
.aggregateFeatures <- function(sparse, featureids) {
    rownames(sparse) <- featureids
    sparse <- sparse[!is.na(rownames(sparse)), , drop = FALSE]
    aggregate.Matrix(sparse, groupings = rownames(sparse), fun = "sum")
}



# Methods ====
#' @rdname aggregateFeatures
#' @export
setMethod(
    "aggregateFeatures",
    signature("bcbioSingleCellANY"),
    function(object) {
        stop("Draft function", call. = FALSE)
        # Reslot the counts into assay and then update the object
    })
