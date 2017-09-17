#' Aggregate Replicates
#'
#' @rdname aggregateReplicates
#' @name aggregateReplicates
#' @family Data Management Utilities
#' @author Rory Kirchner, Michael Steinbaugh
#'
#' @inheritParams AllGenerics
#'
#' @return [bcbioSCDataSet].
NULL



# Constructors ====
.aggregateReplicates <- function(sparse, cellids) {
    tsparse <- t(sparse)
    rownames(tsparse) <- cellids
    aggregate.Matrix(tsparse, cellids, fun = "sum") %>% t
}



# Methods ====
#' @rdname aggregateReplicates
#' @export
setMethod("aggregateReplicates", "bcbioSCDataSet", function(object) {
    stop("Draft function", call. = FALSE)
})
