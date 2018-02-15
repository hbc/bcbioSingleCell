#' Subset Per Sample
#'
#' @rdname subsetPerSample
#' @name subsetPerSample
#' @family Data Management Utilities
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#'
#' @param minCells Minimum number of cells required per sample.
#' @param dir Output directory where to save the subset data.
#'
#' @return Character vector of saved [bcbioSingleCell] subsets.
#'
#' @examples
#' load(system.file(
#'     file.path("extdata", "filtered.rda"),
#'     package = "bcbioSingleCell"))
#'
#' # This will save the subsets per sample to disk
#' subsetPerSample(filtered, dir = "subsetPerSample")
#' dir("subsetPerSample")
#'
#' # Cleanup
#' unlink("subsetPerSample", recursive = TRUE)
NULL



# Methods ======================================================================
#' @rdname subsetPerSample
#' @importFrom dplyr pull
#' @importFrom pbapply pblapply
#' @export
setMethod(
    "subsetPerSample",
    signature("bcbioSingleCell"),
    function(
        object,
        minCells = 200L,
        dir = getwd()) {
        dir.create(dir, recursive = TRUE, showWarnings = FALSE)
        sampleIDs <- sampleMetadata(object) %>%
            pull("sampleID") %>%
            as.character()
        pblapply(seq_along(sampleIDs), function(a) {
            sampleID <- sampleIDs[[a]]
            subset <- selectSamples(
                object,
                sampleID = sampleID)
            # Skip if subset doesn't have enough cells
            if (dim(subset)[[2L]] < minCells) {
                warn(paste(sampleID, "didn't pass minimum cell cutoff"))
                return(NULL)
            }
            assign(sampleID, subset)
            save(
                list = sampleID,
                file = file.path(dir, paste0(sampleID, ".rda")))
            sampleID
        }) %>%
            na.omit() %>%
            as.character()
    })
