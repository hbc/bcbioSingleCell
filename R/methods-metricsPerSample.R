#' Metrics per Sample
#'
#' @rdname metricsPerSample
#' @name metricsPerSample
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#'
#' @return [tibble] containing summary statistics.
#'
#' @examples
#' load(system.file(
#'     file.path("extdata", "bcb.rda"),
#'     package = "bcbioSingleCell"))
#' load(system.file(
#'     file.path("extdata", "seurat.rda"),
#'     package = "bcbioSingleCell"))
#'
#' # bcbioSingleCell
#' metricsPerSample(bcb)
#'
#' # seurat
#' metricsPerSample(seurat)
NULL



# Constructors =================================================================
#' @importFrom dplyr group_by summarize_all
.metricsPerSample <- function(object) {
    metrics(object) %>%
        .[, c(
            "sampleName",
            "nUMI",
            "nGene",
            "nCoding"
        ), drop = FALSE] %>%
        group_by(!!sym("sampleName")) %>%
        summarize_all(sum)
}



# Methods ======================================================================
#' @rdname metricsPerSample
#' @export
setMethod(
    "metricsPerSample",
    signature("bcbioSingleCell"),
    .metricsPerSample)



#' @rdname metricsPerSample
#' @export
setMethod(
    "metricsPerSample",
    signature("seurat"),
    .metricsPerSample)
