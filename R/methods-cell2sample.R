#' Cell to Sample Mappings
#'
#' @name cell2sample
#' @family Data Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#'
#' @examples
#' # bcbioSingleCell ====
#' cell2sample(bcb_small) %>% table()
#'
#' # SingleCellExperiment ====
#' cell2sample(cellranger_small) %>% table()
#'
#' # seurat ====
#' cell2sample(pbmc_small) %>% table()
NULL



# Methods ======================================================================
#' @rdname cell2sample
#' @export
setMethod(
    "cell2sample",
    signature("SingleCellExperiment"),
    function(object) {
        validObject(object)
        cell2sample <- metadata(object)[["cell2sample"]]
        if (!is.factor(cell2sample)) {
            cells <- colnames(object)
            samples <- rownames(sampleData(object))
            cell2sample <- mapCellsToSamples(cells = cells, samples = samples)
        }
        cell2sample %>%
            .[colnames(object)] %>%
            droplevels()
    }
)



#' @rdname cell2sample
#' @export
setMethod(
    "cell2sample",
    signature("seurat"),
    getMethod("cell2sample", "SingleCellExperiment")
)
