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
#' @param envir Environment where to assign the subsets.
#' @param dir Output directory.
#'
#' @return Named vector of saved [bcbioSingleCell] subset file paths.
#'
#' @examples
#' load(system.file("extdata/filtered.rda", package = "bcbioSingleCell"))
#'
#' # This will save the subsets per sample to disk
#' subsetPerSample(filtered, dir = "subsetPerSample")
#' dir_ls("subsetPerSample")
#'
#' # Clean up
#' dir_delete("subsetPerSample")
NULL



# Methods ======================================================================
#' @rdname subsetPerSample
#' @importFrom basejump assignAndSaveData initializeDirectory
#' @importFrom dplyr pull
#' @importFrom fs path_real
#' @importFrom pbapply pblapply
#' @export
setMethod(
    "subsetPerSample",
    signature("bcbioSingleCell"),
    function(
        object,
        minCells = 200L,
        envir = parent.frame(),
        dir = ".") {
        assertIsAnImplicitInteger(minCells)
        assert_all_are_positive(minCells)
        assert_is_environment(envir)
        dir <- initializeDirectory(dir)
        sampleIDs <- sampleData(object) %>%
            pull("sampleID") %>%
            as.character()
        files <- mapply(
            FUN = function(object, sampleID, minCells, envir, dir) {
                subset <- selectSamples(object, sampleID = sampleID)
                # Skip if subset doesn't have enough cells
                if (ncol(subset) < minCells) {
                    warn(paste(sampleID, "didn't pass minimum cell cutoff"))
                    return(NULL)
                }
                assignAndSaveData(
                    name = sampleID,
                    object = subset,
                    dir = dir)
            },
            sampleID = sampleIDs,
            MoreArgs = list(
                object = object,
                minCells = minCells,
                envir = envir,
                dir = dir
            ),
            SIMPLIFY = FALSE,
            USE.NAMES = TRUE)
        # Drop any `NULL` subsets
        files <- Filter(Negate(is.null), files)
        names <- names(files)
        files <- unlist(files)
        files <- path_real(files)
        names(files) <- names
        files
    })
