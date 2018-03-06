#' Cell to Sample Mappings
#'
#' @rdname cell2sample
#' @name cell2sample
#'
#' @inheritParams general
#'
#' @examples
#' # bcbioSingleCell
#' cell2sample(bcb) %>% table()
#'
#' # seurat
#' cell2sample(pbmc_small) %>% table()
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
            abort(paste(
                "cell2sample is not a factor.",
                "Run `updateObject()`."
            ))
        }
        cell2sample %>%
            .[colnames(object)] %>%
            droplevels()
    })



#' @rdname cell2sample
#' @export
setMethod(
    "cell2sample",
    signature("seurat"),
    function(object) {
        # Attempt to use stashed factor first
        cell2sample <- bcbio(object, "cell2sample")
        if (!is.factor(cell2sample)) {
            cells <- colnames(slot(object, "data"))
            samples <- rownames(sampleData(object))
            cell2sample <- mapCellsToSamples(cells = cells, samples = samples)
        }
        cell2sample
    })
