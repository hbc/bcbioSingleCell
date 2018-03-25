#' Cell to Sample Mappings
#'
#' @name cell2sample
#' @family Data Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#'
#' @examples
#' load(system.file("extdata/bcb_small.rda", package = "bcbioSingleCell"))
#'
#' # bcbioSingleCell ====
#' cell2sample(bcb_small) %>% table()
#'
#' # seurat ====
#' cell2sample(pbmc_small) %>% table()
NULL



# Constructors =================================================================
.cell2sample <- function(object) {
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
