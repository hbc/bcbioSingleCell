#' Cell to Sample Mappings
#'
#' @name cell2sample
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#'
#' @examples
#' # bcbioSingleCell ====
#' cell2sample(bcb_small) %>% table()
#'
#' # seurat ====
#' cell2sample(pbmc_small) %>% table()
#' cell2sample(seurat_small) %>% table()
NULL



# Constructors =================================================================
.cell2sample <- function(object) {
    validObject(object)
    cell2sample <- metadata(object)[["cell2sample"]]
    if (!is.factor(cell2sample)) {
        warn(paste(
            "cell2sample is not stashed in object.",
            "Calculating on the fly instead."
        ))
        # FIXME Add support for `colnames()` to seurat
        cells <- rownames(colData(object))
        samples <- rownames(sampleData(object))
        cell2sample <- mapCellsToSamples(cells = cells, samples = samples)
    }
    cell2sample %>%
        # FIXME Add support for `colnames()` to seurat
        .[rownames(colData(object))] %>%
        droplevels()
}



# Methods ======================================================================
#' @rdname cell2sample
#' @export
setMethod(
    "cell2sample",
    signature("bcbioSingleCell"),
    .cell2sample
)



#' @rdname cell2sample
#' @export
setMethod(
    "cell2sample",
    signature("seurat"),
    .cell2sample
)
