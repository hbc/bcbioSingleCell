#' Cell to Sample Mappings
#'
#' @rdname cell2sample
#' @name cell2sample
#' @author Michael Steinbaugh
#'
#' @inheritParams AllGenerics
NULL



# Methods ====
#' @rdname cell2sample
#' @export
setMethod(
    "cell2sample",
    signature("bcbioSingleCell"),
    function(object) {
        # Define the cell2sample mappings
        # This uses a stashed `data.frame` as of v0.0.22, for better speed
        cell2sample <- metadata(object)[["cell2sample"]]
        if (is.null(cell2sample)) {
            cell2sample <- .cell2sample(
                cells = colData[["cellID"]],
                samples = metadata[["sampleID"]]
            )
        }
        cell2sample[["sampleID"]] <- as.factor(cell2sample[["sampleID"]])
        cell2sample
    })
