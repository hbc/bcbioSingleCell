#' Aggregate Features
#'
#' @rdname aggregateFeatures
#' @name aggregateFeatures
#' @family Data Management Utilities
#' @author Rory Kirchner, Michael Steinbaugh
#'
#' @importFrom basejump aggregateFeatures
#'
#' @inherit basejump::aggregateFeatures
NULL



# Methods ====
#' @rdname aggregateFeatures
#' @export
setMethod(
    "aggregateFeatures",
    signature("bcbioSingleCell"),
    function(object) {
        warning(paste(
            "Draft function.",
            "Returning an aggregated counts matrix."
        ), call. = FALSE)
        aggregateFeatures(assay(object), features = features)
    })
