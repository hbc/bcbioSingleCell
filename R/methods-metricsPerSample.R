#' Metrics per Sample
#'
#' @name metricsPerSample
#' @family Quality Control Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#'
#' @return `tbl_df` containing summary statistics.
#'
#' @examples
#' load(system.file("extdata/bcb_small.rda", package = "bcbioSingleCell"))
#' load(system.file("extdata/seurat_small.rda", package = "bcbioSingleCell"))
#'
#' # bcbioSingleCell ====
#' metricsPerSample(bcb_small)
#'
#' # seurat ====
#' metricsPerSample(seurat_small)
NULL



# Constructors =================================================================
.metricsPerSample <- function(object) {
    metrics(object) %>%
        .[, c("sampleName", "nUMI", "nGene", "nCoding"), drop = FALSE] %>%
        group_by(!!sym("sampleName")) %>%
        summarize_all(sum)
}



# Methods ======================================================================
#' @rdname metricsPerSample
#' @export
setMethod(
    "metricsPerSample",
    signature("bcbioSingleCell"),
    .metricsPerSample
)



#' @rdname metricsPerSample
#' @export
setMethod(
    "metricsPerSample",
    signature("seurat"),
    .metricsPerSample
)
