#' Cell to Sample Mappings
#'
#' @name cell2sample
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#'
#' @examples
#' # bcbioSingleCell
#' cell2sample(bcb_small) %>% table()
#'
#' # seurat
#' cell2sample(pbmc_small) %>% table()
#' cell2sample(seurat_small) %>% table()
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
        cell2sample <- metadata(object)[["cell2sample"]]
        if (!is.factor(cell2sample)) {
            cells <- colnames(slot(object, "data"))
            samples <- rownames(sampleData(object))
            cell2sample <- mapCellsToSamples(cells = cells, samples = samples)
        }
        assert_is_subset(names(cell2sample), rownames(colData(object)))
        cell2sample %>%
            .[rownames(colData(object))] %>%
            droplevels()
    }
)
