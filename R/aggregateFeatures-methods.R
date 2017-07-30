#' Aggregate Features
#'
#' @rdname aggregateReplicates
#' @name aggregateReplicates
#' @author Rory Kirchner, Michael Steinbaugh
#'
#' @return [bcbioSCDataSet].
NULL



# Constructors ====
.aggregateFeatures <- function(sparse, featureids) {
    rownames(sparse) <- featureids
    sparse <- sparse[!is.na(rownames(sparse)), ]
    aggregate.Matrix(sparse, rownames(sparse), fun = "sum")
}



# Methods ====
#' @rdname aggregateFeatures
#' @export
setMethod("aggregateFeatures", "bcbioSCDataSet", function(object) {
    stop("Draft function", call. = FALSE)
})
