#' Aggregate Replicates
#'
#' @rdname aggregateReplicates
#' @name aggregateReplicates
#' @family Data Management Utilities
#' @author Rory Kirchner, Michael Steinbaugh
#' @keywords internal
#'
#' @inheritParams AllGenerics
#'
#' @param object Sparse counts matrix (e.g. `dgCMatrix`).
#' @param cellids Cellular barcode identifiers.
#'
#' @return `dgCMatrix`.
NULL



# Constructors ====
#' @importFrom Matrix.utils aggregate.Matrix
.aggregateSparseReplicates <- function(object, cellids) {
    tsparse <- t(object)
    rownames(tsparse) <- cellids
    tsparse %>%
        aggregate.Matrix(groupings = cellids, fun = "sum") %>%
        t()
}



# Methods ====
#' @rdname aggregateReplicates
#' @export
setMethod(
    "aggregateReplicates",
    signature("bcbioSingleCell"),
    function(object, cellids) {
        warning(paste(
            "Draft function.",
            "Returning an aggregated counts matrix."
        ), call. = FALSE)
        .aggregateSparseReplicates(
            object = assay(object),
            cellids = cellids
        )
    })



#' @rdname aggregateReplicates
#' @export
setMethod(
    "aggregateReplicatess",
    signature("dgCMatrix"),
    .aggregateSparseReplicates)
