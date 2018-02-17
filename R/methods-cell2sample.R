#' Cell to Sample Mappings
#'
#' @rdname cell2sample
#' @name cell2sample
#'
#' @inheritParams general
NULL



# Methods ======================================================================
#' @rdname cell2sample
#' @export
setMethod(
    "cell2sample",
    signature("bcbioSingleCell"),
    function(object) {
        cell2sample <- metadata(object)[["cell2sample"]]

        if (!is.factor(cell2sample)) {
            cell2sample <- mapCellsToSamples(
                cells = colnames(object),
                samples = rownames(sampleMetadata(object))
            )
            return(cell2sample)
        }

        # Version-specific fixes
        if (metadata(object)[["version"]] == "0.0.22") {
            # v0.0.22 stashed this as a data.frame instead
            if (is.data.frame(cell2sample)) {
                cells <- as.character(cell2sample[["cellID"]])
                samples <- as.factor(cell2sample[["sampleID"]])
                cell2sample <- samples
                names(cell2sample) <- cells
            } else {
                cell2sample <- NULL
            }
        }

        cell2sample <- cell2sample[colnames(object)]
        cell2sample <- droplevels(cell2sample)
        cell2sample
    })



#' @rdname cell2sample
#' @export
setMethod(
    "cell2sample",
    signature("seurat"),
    function(object) {
        cell2sample <- bcbio(object, "cell2sample")
        if (!is.factor(cell2sample)) {
            cell2sample <- mapCellsToSamples(
                cells = colnames(slot(object, "data")),
                samples = rownames(sampleMetadata(object)))
        }
        cell2sample
    })
