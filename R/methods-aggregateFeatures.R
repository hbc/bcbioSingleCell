# FIXME This function needs to be updated



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
#' @inheritParams AllGenerics
NULL



# Methods ====
#' @rdname aggregateFeatures
#' @export
setMethod(
    "aggregateFeatures",
    signature("bcbioSingleCell"),
    function(object, features) {
        warning(paste(
            "Draft function.",
            "Returning a only the aggregated counts matrix."
        ), call. = FALSE)
        aggregateFeatures(assay(object), features = features)
    })
