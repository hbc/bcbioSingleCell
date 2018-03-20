#' Metrics per Sample
#'
#' @name metricsPerSample
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#'
#' @return `tbl_df` containing summary statistics.
#'
#' @examples
#' # bcbioSingleCell ====
#' metricsPerSample(bcb_small)
#'
#' # seurat ====
#' metricsPerSample(seurat_small)
NULL



# Constructors =================================================================
#' @importFrom dplyr group_by summarize_all
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
