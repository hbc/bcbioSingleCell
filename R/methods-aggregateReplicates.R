#' Aggregate Replicates
#'
#' @rdname aggregateReplicates
#' @name aggregateReplicates
#' @family Data Management Utilities
#' @author Rory Kirchner, Michael Steinbaugh
#'
#' @inheritParams AllGenerics
#'
#' @return [bcbioSingleCell].
NULL



# Constructors ====
#' Aggregate Replicates Constructor
#'
#' @author Rory Kirchner
#' @keywords internal
#' @noRd
#'
#' @param sparse Sparse counts matrix (e.g. `dgCMatrix`).
#' @param cellids Cellular barcode identifiers.
#'
#' @return [dgCMatrix].
.aggregateReplicates <- function(sparse, cellids) {
    tsparse <- t(sparse)
    rownames(tsparse) <- cellids
    tsparse %>%
        aggregate.Matrix(groupings = cellids, fun = "sum") %>%
        t()
}



# Methods ====
#' @rdname aggregateReplicates
#' @export
setMethod("aggregateReplicates", "bcbioSingleCell", function(object) {
    stop("Draft function", call. = FALSE)
    # Reslot the counts into assay and then update the object.
    # Slot `sampleNameAggregate` into `sampleMetadata()`.
})
